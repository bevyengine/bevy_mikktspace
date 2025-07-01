#![cfg_attr(docsrs, feature(doc_auto_cfg))]
#![doc(
    html_logo_url = "https://bevy.org/assets/icon.png",
    html_favicon_url = "https://bevy.org/assets/icon.png"
)]
#![no_std]

//! Mikktspace ported to *safe* Rust.
//!
//! This library is based on Morten S. Mikkelsen's original tangent space algorithm
//! implementation, written in C. The original source code can be found at
//! <https://archive.blender.org/wiki/index.php/Dev:Shading/Tangent_Space_Normal_Maps>
//! and includes the following licence:
//!
//! > Copyright (C) 2011 by Morten S. Mikkelsen
//! >
//! > This software is provided 'as-is', without any express or implied
//! > warranty.  In no event will the authors be held liable for any damages
//! > arising from the use of this software.
//! >
//! > Permission is granted to anyone to use this software for any purpose,
//! > including commercial applications, and to alter it and redistribute it
//! > freely, subject to the following restrictions:
//! >
//! > 1. The origin of this software must not be misrepresented; you must not
//! > claim that you wrote the original software. If you use this software
//! > in a product, an acknowledgment in the product documentation would be
//! > appreciated but is not required.
//! >
//! > 2. Altered source versions must be plainly marked as such, and must not be
//! > misrepresented as being the original software.
//! >
//! > 3. This notice may not be removed or altered from any source distribution.

extern crate alloc;

use alloc::{vec, vec::Vec};
use bitflags::bitflags;

mod geometry;
mod math;

#[cfg(test)]
mod tests;

use geometry::GeometryExt;

pub use geometry::Geometry;
pub use math::VectorSpace;

use crate::math::{abs, not_zero};

/// Generates tangents for the input geometry.
pub fn generate_tangents<G: Geometry>(geometry: &mut G) {
    // This allows any invocation of `generate_tangents` to be compared against
    // the reference C implementation without any additional tooling.
    #[cfg(test)]
    let geometry = &mut tests::shadow::Shadow::new(geometry);

    generate_tangents_with(geometry, -1.0)
}

/// Generates tangents for the input geometry and a given `cos(angular_threshold)`.
fn generate_tangents_with<G: Geometry>(geometry: &mut G, thres_cos: f32) {
    let (mut triangles, mut vertices, tspaces_count) = generate_initial_vertices(geometry);

    weld_vertices(geometry, &mut vertices);

    let info = split_degenerate(geometry, &mut triangles, &mut vertices);

    init_tri_info(geometry, info.good_triangles, info.good_vertices);

    // Groups allocate their own buffers out of this arena
    let mut group_buffer = vec![0; info.good_vertices.len()];

    let mut groups =
        build_4_rule_groups(info.good_triangles, &mut group_buffer, info.good_vertices);

    let mut tspaces = generate_tspaces(
        geometry,
        tspaces_count,
        info.good_triangles,
        &mut groups,
        info.good_vertices,
        thres_cos,
    );

    approximate_degenerate_triangles(geometry, &mut tspaces, info);

    let mut tspaces = tspaces.into_iter();
    for face in 0..geometry.num_faces() {
        let verts = geometry.num_vertices_of_face(face);
        if verts == 3 || verts == 4 {
            // I've decided to let degenerate triangles and group-with-anythings
            // vary between left/right hand coordinate systems at the vertices.
            // All healthy triangles on the other hand are built to always be either or.

            // set data
            for (vert, tspace) in (0..verts).zip(&mut tspaces) {
                geometry.set_tangent(
                    tspace.os,
                    tspace.ot,
                    tspace.mag_s,
                    tspace.mag_t,
                    tspace.orient,
                    face,
                    vert,
                );
            }
        }
    }
}

struct TriangleInfo<V: VectorSpace> {
    face_neighbors: [Option<usize>; 3],
    assigned_group: [Option<usize>; 3],
    os: V::Vec3,
    ot: V::Vec3,
    mag_s: f32,
    mag_t: f32,
    /// Index of the face this triangle maps to, in the original faces.
    original_face_index: usize,
    flags: TriangleFlags,
    /// Offset of the first vertex of this triangle, in the original vertices.
    vertex_offset: usize,
    /// Offsets of the vertices of this triangle, relative to the triangle index
    vertex_indices: [u8; 3],
}

impl<V: VectorSpace> Clone for TriangleInfo<V>
where
    V::Vec3: Copy,
{
    fn clone(&self) -> Self {
        *self
    }
}

impl<V: VectorSpace> Copy for TriangleInfo<V> where V::Vec3: Copy {}

bitflags! {
    #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
    struct TriangleFlags: u32 {
        const MARK_DEGENERATE = 0b00000001;
        const QUAD_ONE_DEGEN_TRI = 0b00000010;
        const GROUP_WITH_ANY = 0b00000100;
        const ORIENT_PRESERVING = 0b00001000;
    }
}

impl<V: VectorSpace> TriangleInfo<V> {
    fn zero() -> Self {
        Self {
            face_neighbors: [None; 3],
            assigned_group: [None; 3],
            os: [0.; 3].into(),
            ot: [0.; 3].into(),
            mag_s: 0.0,
            mag_t: 0.0,
            original_face_index: 0,
            flags: TriangleFlags::empty(),
            vertex_offset: 0,
            vertex_indices: [0; 3],
        }
    }
}

/// Generate initial triangles and vertex indices, from original geometry.
fn generate_initial_vertices<G: Geometry>(
    geometry: &G,
) -> (Vec<TriangleInfo<G::Space>>, Vec<usize>, usize) {
    let triangles_count = geometry
        .face_indices()
        .map(|face| geometry.num_vertices_of_face(face))
        .map(|verts| match verts {
            3 => 1,
            4 => 2,
            // Unsupported face-types
            _ => 0,
        })
        .sum::<usize>();

    let mut triangles_info = Vec::with_capacity(triangles_count);
    let mut vertex_indices = Vec::with_capacity(3 * triangles_count);

    let mut vertex_offset = 0;

    for original_face_index in geometry.face_indices() {
        let mut face_vertex_indices = geometry.vertex_indices(original_face_index);

        let mut tri_info_a = TriangleInfo {
            original_face_index,
            vertex_offset,
            ..TriangleInfo::zero()
        };

        vertex_offset += face_vertex_indices.len();

        match face_vertex_indices.len() {
            3 => {
                // For tris
                tri_info_a.vertex_indices = [0, 1, 2];

                triangles_info.push(tri_info_a);
                vertex_indices.extend(face_vertex_indices);
            }
            4 => {
                // For quads
                let mut tri_info_b = tri_info_a;

                let i = [(); 4].map(|_| face_vertex_indices.next().unwrap());
                let t = i.map(|i| geometry.tex_coord_by_index(i));

                // Figure out the best cut for the quad
                let quad_diagonal_is_02 = match G::Space::distance_squared(t[2], t[0])
                    .total_cmp(&G::Space::distance_squared(t[3], t[1]))
                {
                    core::cmp::Ordering::Less => true,
                    core::cmp::Ordering::Greater => false,
                    core::cmp::Ordering::Equal => {
                        let p = i.map(|i| geometry.position_by_index(i));
                        G::Space::distance_squared(p[3], p[1])
                            >= G::Space::distance_squared(p[2], p[0])
                    }
                };

                // Apply indices for the cut we determined
                if quad_diagonal_is_02 {
                    tri_info_a.vertex_indices = [0, 1, 2];
                    tri_info_b.vertex_indices = [0, 2, 3];
                } else {
                    tri_info_a.vertex_indices = [0, 1, 3];
                    tri_info_b.vertex_indices = [1, 2, 3];
                }

                triangles_info.push(tri_info_a);
                triangles_info.push(tri_info_b);
                vertex_indices.extend_from_slice(&tri_info_a.vertex_indices.map(|i| i as usize));
                vertex_indices.extend_from_slice(&tri_info_b.vertex_indices.map(|i| i as usize));
            }
            _ => {
                // Only generate for tris or quads
                continue;
            }
        }
    }

    (triangles_info, vertex_indices, vertex_offset)
}

struct TSpace<V: VectorSpace> {
    os: V::Vec3,
    mag_s: f32,
    ot: V::Vec3,
    mag_t: f32,
    counter: u8,
    orient: bool,
}

impl<V: VectorSpace> Clone for TSpace<V>
where
    V::Vec3: Copy,
{
    fn clone(&self) -> Self {
        *self
    }
}

impl<V: VectorSpace> Copy for TSpace<V> where V::Vec3: Copy {}

impl<V: VectorSpace> TSpace<V> {
    fn zero() -> Self {
        Self {
            os: [0.; 3].into(),
            mag_s: 0.0,
            ot: [0.; 3].into(),
            mag_t: 0.0,
            counter: 0,
            orient: false,
        }
    }

    fn mean(self, other: Self) -> Self {
        if self.mag_s == other.mag_s
            && self.mag_t == other.mag_t
            && self.os == other.os
            && self.ot == other.ot
        {
            Self {
                mag_s: self.mag_s,
                mag_t: self.mag_t,
                os: self.os,
                ot: self.ot,
                ..Self::zero()
            }
        } else {
            Self {
                mag_s: 0.5 * (self.mag_s + other.mag_s),
                mag_t: 0.5 * (self.mag_t + other.mag_t),
                os: V::normalize_or_zero(self.os + other.os),
                ot: V::normalize_or_zero(self.ot + other.ot),
                ..Self::zero()
            }
        }
    }
}

// To avoid visual errors (distortions/unwanted hard edges in lighting), when using sampled normal
// maps, the normal map sampler must use the exact inverse of the pixel shader transformation.
// The most efficient transformation we can possibly do in the pixel shader is
// achieved by using, directly, the "unnormalized" interpolated tangent, bitangent and vertex
// normal: vT, vB and vN.
// pixel shader (fast transform out)
// vNout = normalize( vNt.x * vT + vNt.y * vB + vNt.z * vN );
// where vNt is the tangent space normal. The normal map sampler must likewise use the
// interpolated and "unnormalized" tangent, bitangent and vertex normal to be compliant with the
// pixel shader.
// sampler does (exact inverse of pixel shader):
// float3 row0 = cross(vB, vN);
// float3 row1 = cross(vN, vT);
// float3 row2 = cross(vT, vB);
// float fSign = dot(vT, row0)<0 ? -1 : 1;
// vNt = normalize( fSign * float3(dot(vNout,row0), dot(vNout,row1), dot(vNout,row2)) );
// where vNout is the sampled normal in some chosen 3D space.
//
// Should you choose to reconstruct the bitangent in the pixel shader instead
// of the vertex shader, as explained earlier, then be sure to do this in the normal map sampler
// also.
// Finally, beware of quad triangulations. If the normal map sampler doesn't use the same
// triangulation of
// quads as your renderer then problems will occur since the interpolated tangent spaces will differ
// eventhough the vertex level tangent spaces match. This can be solved either by triangulating
// before
// sampling/exporting or by using the order-independent choice of diagonal for splitting quads
// suggested earlier.
// However, this must be used both by the sampler and your tools/rendering pipeline.
// internal structure

#[derive(Copy, Clone)]
struct Group<'a> {
    face_indices: &'a [usize],
    vertex_representative: usize,
    orient_preserving: bool,
}

#[derive(Clone, PartialEq)]
struct SubGroup {
    faces_count: usize,
    tri_members: Vec<usize>,
}

impl SubGroup {
    const fn zero() -> Self {
        Self {
            faces_count: 0,
            tri_members: Vec::new(),
        }
    }
}

#[derive(Copy, Clone, PartialEq, Eq)]
struct Edge {
    /// Index of the starting vertex
    i0: usize,
    /// Index of the ending vertex
    i1: usize,
    /// Index of the face this edge belongs to
    face: usize,
}

impl Ord for Edge {
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        use core::cmp::Ordering::Equal;

        match self.i0.cmp(&other.i0) {
            Equal => match self.i1.cmp(&other.i1) {
                Equal => self.face.cmp(&other.face),
                ord => ord,
            },
            ord => ord,
        }
    }
}

impl PartialOrd for Edge {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// Set the tangent space values of degenerate triangles based on connected good triangles.
///
/// Degenerate quads with one good triangle will be fixed by copying a space from
/// the good triangle to the coinciding vertex.
/// All other degenerate triangles will just copy a space from any good triangle
/// with the same welded index.
fn approximate_degenerate_triangles<G: Geometry>(
    geometry: &G,
    tspaces: &mut [TSpace<G::Space>],
    info: Info<'_, G::Space>,
) {
    // deal with degenerate triangles
    // punishment for degenerate triangles is O(N^2)
    // degenerate triangles on a quad with one good triangle are skipped
    // here but processed in the next loop
    let full_bad = info
        .degenerate_triangles
        .iter()
        .enumerate()
        .filter(|(_, triangle)| !triangle.flags.contains(TriangleFlags::QUAD_ONE_DEGEN_TRI))
        .flat_map(|(t, a)| (0..3).map(move |i| (t, a, i)))
        .map(|(t, a, i)| (info.degenerate_vertices[t * 3 + i], a, i))
        .filter_map(|(degen_index, a, i)| {
            info.good_vertices
                .iter()
                .enumerate()
                .find(|&(_, &index)| index == degen_index)
                .map(|(j, _)| (i, j, a))
        })
        .map(|(i, j, a)| {
            let b = &info.good_triangles[j / 3];
            let j = j % 3;

            let dst = a.vertex_indices[i] as usize + a.vertex_offset;
            let src = b.vertex_indices[j] as usize + b.vertex_offset;

            (dst, src)
        });

    // deal with degenerate quads with one good triangle
    // this triangle belongs to a quad where the
    // other triangle is degenerate
    let partial_bad = info
        .degenerate_triangles
        .iter()
        .filter(|triangle| triangle.flags.contains(TriangleFlags::QUAD_ONE_DEGEN_TRI))
        .filter_map(|info| {
            let a = (0..4)
                .find(|index| !info.vertex_indices.contains(index))
                .unwrap() as usize;

            info.vertex_indices
                .iter()
                .map(|&b| b as usize)
                .find(|&b| {
                    geometry.position(info.original_face_index, a)
                        == geometry.position(info.original_face_index, b)
                })
                .map(|b| (info.vertex_offset + a, info.vertex_offset + b))
        });

    for (dst, src) in full_bad.chain(partial_bad) {
        tspaces[dst] = tspaces[src];
    }
}

/// Make tspaces, each group is split up into subgroups if necessary
/// based on `angular_threshold`.
/// Finally a tangent space is made for every resulting subgroup.
fn generate_tspaces<G: Geometry>(
    geometry: &G,
    tspaces_count: usize,
    triangles_info: &[TriangleInfo<G::Space>],
    groups: &mut [Group],
    vertex_indices: &[usize],
    thres_cos: f32,
) -> Vec<TSpace<G::Space>> {
    let mut tspaces = vec![
        TSpace {
            os: [1.0, 0.0, 0.0].into(),
            mag_s: 1.0,
            ot: [0.0, 1.0, 0.0].into(),
            mag_t: 1.0,
            ..TSpace::zero()
        };
        tspaces_count
    ];

    let max_faces_count = groups
        .iter()
        .map(|group| group.face_indices.len())
        .max()
        .unwrap_or(0);

    if max_faces_count == 0 {
        return tspaces;
    }

    let mut sub_group_tspace = vec![TSpace::zero(); max_faces_count];
    let mut uni_sub_groups = vec![SubGroup::zero(); max_faces_count];
    let mut tmp_members = vec![0; max_faces_count];

    for (g, group) in groups.iter_mut().enumerate() {
        let mut unique_sub_groups = 0;

        let face_indices = group.face_indices;

        for &f in face_indices {
            let a = &triangles_info[f];

            let index = a
                .assigned_group
                .iter()
                .position(|&group| group == Some(g))
                .unwrap();

            let vertex_index = vertex_indices[f * 3 + index];
            let n = geometry.normal_by_index(vertex_index);
            let a_v_os = G::Space::normalize_or_zero(G::Space::project(a.os, n));
            let a_v_ot = G::Space::normalize_or_zero(G::Space::project(a.ot, n));

            let members = face_indices
                .iter()
                .copied()
                .filter(|&t| {
                    let b = &triangles_info[t];
                    let b_v_os = G::Space::normalize_or_zero(G::Space::project(b.os, n));
                    let b_v_ot = G::Space::normalize_or_zero(G::Space::project(b.ot, n));

                    let flags = a.flags | b.flags;
                    let any = flags.contains(TriangleFlags::GROUP_WITH_ANY);
                    // make sure triangles which belong to the same quad are joined.
                    let same_org_face = a.original_face_index == b.original_face_index;
                    let cos_s = G::Space::dot(a_v_os, b_v_os);
                    let cos_t = G::Space::dot(a_v_ot, b_v_ot);
                    any || same_org_face || cos_s > thres_cos && cos_t > thres_cos
                })
                .fold(0, |members, t| {
                    tmp_members[members] = t;
                    members + 1
                });

            tmp_members[0..members].sort();

            let tmp_group = SubGroup {
                faces_count: members,
                tri_members: tmp_members.clone(),
            };

            let l = uni_sub_groups
                .iter()
                .take(unique_sub_groups)
                .position(|group| group == &tmp_group);

            if l.is_none() {
                uni_sub_groups[unique_sub_groups].faces_count = members;
                uni_sub_groups[unique_sub_groups].tri_members = tmp_group.tri_members.clone();

                sub_group_tspace[unique_sub_groups] = eval_tspace(
                    geometry,
                    &tmp_group.tri_members[..members],
                    vertex_indices,
                    triangles_info,
                    group.vertex_representative,
                );
                unique_sub_groups += 1;
            }

            let sub_tspace = sub_group_tspace[l.unwrap_or(unique_sub_groups - 1)];

            let index = a.vertex_indices[index] as usize + a.vertex_offset;

            let tspaces_out = &mut tspaces[index];
            let sub_tspace = if tspaces_out.counter == 1 {
                tspaces_out.mean(sub_tspace)
            } else {
                sub_tspace
            };

            *tspaces_out = TSpace {
                counter: tspaces_out.counter + 1,
                orient: group.orient_preserving,
                ..sub_tspace
            };
        }
    }

    tspaces
}

fn eval_tspace<G: Geometry>(
    geometry: &G,
    face_indices_buffer: &[usize],
    vertices: &[usize],
    triangles_info: &[TriangleInfo<G::Space>],
    vertex_representative: usize,
) -> TSpace<G::Space> {
    let chunks = vertices.chunks_exact(3).zip(triangles_info);

    let (mut res, angle_sum) = face_indices_buffer
        .iter()
        .map(|&f| chunks.clone().nth(f).unwrap())
        .filter(|(_vertex_indices, info)| {
            // only valid triangles get to add their contribution
            !info.flags.contains(TriangleFlags::GROUP_WITH_ANY)
        })
        .fold(
            (TSpace::zero(), 0.),
            |(mut res, angle_sum), (vertices, info)| {
                let this = vertices
                    .iter()
                    .position(|&v| v == vertex_representative)
                    .unwrap();

                let mut vertices = vertices.iter().cycle().skip(this + 2).copied();
                let i = [(); 3].map(|_| vertices.next().unwrap());

                let n = geometry.normal_by_index(i[1]);

                let angle = {
                    let p = i.map(|i| geometry.position_by_index(i));
                    let l01 = G::Space::project(p[0] - p[1], n);
                    let l21 = G::Space::project(p[2] - p[1], n);
                    G::Space::angle_between(l01, l21)
                };

                let v_os = G::Space::normalize_or_zero(G::Space::project(info.os, n));
                let v_ot = G::Space::normalize_or_zero(G::Space::project(info.ot, n));

                res.os = res.os + G::Space::scale(v_os, angle);
                res.ot = res.ot + G::Space::scale(v_ot, angle);
                res.mag_s += angle * info.mag_s;
                res.mag_t += angle * info.mag_t;

                (res, angle_sum + angle)
            },
        );

    res.os = G::Space::normalize_or_zero(res.os);
    res.ot = G::Space::normalize_or_zero(res.ot);

    if angle_sum > 0. {
        res.mag_s /= angle_sum;
        res.mag_t /= angle_sum;
    }

    res
}

/// Based on the 4 rules, identify groups based on connectivity.
fn build_4_rule_groups<'a, V: VectorSpace>(
    triangles_info: &mut [TriangleInfo<V>],
    mut face_indices_buffer: &'a mut [usize],
    vertices: &[usize],
) -> Vec<Group<'a>> {
    let mut group_index = 0;
    let mut groups = Vec::with_capacity(vertices.len());

    for (index, &vertex_representative) in vertices.iter().enumerate() {
        let f = index / 3;
        let i = index % 3;

        let info = &mut triangles_info[f];

        if info.flags.contains(TriangleFlags::GROUP_WITH_ANY) || info.assigned_group[i].is_some() {
            continue;
        }

        info.assigned_group[i] = Some(group_index);
        let orient_preserving = info.flags.contains(TriangleFlags::ORIENT_PRESERVING);

        face_indices_buffer[0] = f;
        let mut face_indices_len = 1;

        let or_pre = info.flags.contains(TriangleFlags::ORIENT_PRESERVING);

        for neighbor in [i, i + 2]
            .map(|i| info.face_neighbors[i % 3])
            .into_iter()
            .flatten()
        {
            let answer = assign_recur(
                vertices,
                triangles_info,
                neighbor,
                group_index,
                face_indices_buffer,
                vertex_representative,
                orient_preserving,
                &mut face_indices_len,
            );
            let or_pre2 = triangles_info[neighbor]
                .flags
                .contains(TriangleFlags::ORIENT_PRESERVING);
            let diff = or_pre != or_pre2;
            assert!(answer || diff);
        }

        group_index += 1;

        let face_indices;
        (face_indices, face_indices_buffer) = face_indices_buffer.split_at_mut(face_indices_len);

        groups.push(Group {
            face_indices,
            vertex_representative,
            orient_preserving,
        });
    }

    groups
}

#[expect(clippy::too_many_arguments)]
fn assign_recur<V: VectorSpace>(
    vertex_indices: &[usize],
    triangles_info: &mut [TriangleInfo<V>],
    my_tri_index: usize,
    group_index: usize,
    face_indices_buffer: &mut [usize],
    vertex_representative: usize,
    orient_preserving: bool,
    face_indices_len: &mut usize,
) -> bool {
    let my_tri_info = &mut triangles_info[my_tri_index];

    // track down vertex
    let i = vertex_indices
        .iter()
        .copied()
        .skip(3 * my_tri_index)
        .take(3)
        .position(|v| v == vertex_representative)
        .unwrap();

    // early out
    if let Some(existing_group) = my_tri_info.assigned_group[i] {
        return group_index == existing_group;
    }

    // first to group with a group-with-anything triangle
    // determines it's orientation.
    // This is the only existing order dependency in the code!!
    if my_tri_info.flags.contains(TriangleFlags::GROUP_WITH_ANY)
        && my_tri_info.assigned_group.iter().all(|g| g.is_none())
    {
        my_tri_info
            .flags
            .set(TriangleFlags::ORIENT_PRESERVING, orient_preserving);
    }

    if my_tri_info.flags.contains(TriangleFlags::ORIENT_PRESERVING) != orient_preserving {
        return false;
    }

    face_indices_buffer[*face_indices_len] = my_tri_index;
    *face_indices_len += 1;

    my_tri_info.assigned_group[i] = Some(group_index);

    for neighbor in [i, i + 2]
        .map(|i| my_tri_info.face_neighbors[i % 3])
        .into_iter()
        .flatten()
    {
        assign_recur(
            vertex_indices,
            triangles_info,
            neighbor,
            group_index,
            face_indices_buffer,
            vertex_representative,
            orient_preserving,
            face_indices_len,
        );
    }

    true
}

/// Evaluate triangle level attributes and neighbor list.
fn init_tri_info<G: Geometry>(
    geometry: &G,
    triangles: &mut [TriangleInfo<G::Space>],
    vertices: &[usize],
) {
    // evaluate first order derivatives
    triangles
        .iter_mut()
        .map(|info| {
            *info = TriangleInfo {
                original_face_index: info.original_face_index,
                vertex_offset: info.vertex_offset,
                vertex_indices: info.vertex_indices,
                flags: info.flags | TriangleFlags::GROUP_WITH_ANY,
                ..TriangleInfo::zero()
            };

            info
        })
        .zip(vertices.chunks_exact(3))
        .filter_map(|(info, i)| {
            let i = [i[0], i[1], i[2]];
            let v = i.map(|i| geometry.position_by_index(i));
            let t = i.map(|i| geometry.tex_coord_by_index(i).into());
            let a = t[1][0] - t[0][0];
            let b = t[1][1] - t[0][1];
            let c = t[2][0] - t[0][0];
            let d = t[2][1] - t[0][1];
            let d1 = v[1] - v[0];
            let d2 = v[2] - v[0];

            let signed_area = a * d - b * c;
            let os = G::Space::scale(d1, d) - G::Space::scale(d2, b);
            let ot = G::Space::scale(d1, -c) + G::Space::scale(d2, a);

            info.flags.set(
                TriangleFlags::ORIENT_PRESERVING,
                signed_area.is_sign_positive(),
            );

            not_zero(signed_area).then_some((info, signed_area, os, ot))
        })
        .map(|(info, signed_area_stx2, os, ot)| {
            info.os = G::Space::normalize_or_zero(os);
            info.ot = G::Space::normalize_or_zero(ot);

            if !info.flags.contains(TriangleFlags::ORIENT_PRESERVING) {
                info.os = -info.os;
                info.ot = -info.ot;
            }

            // evaluate magnitudes prior to normalization of vOs and vOt
            let abs_area = abs(signed_area_stx2);
            info.mag_s = G::Space::length(os) / abs_area;
            info.mag_t = G::Space::length(ot) / abs_area;

            info
        })
        .filter(|info| not_zero(info.mag_s) && not_zero(info.mag_t))
        .for_each(|info| info.flags.remove(TriangleFlags::GROUP_WITH_ANY));

    triangles
        .chunk_by_mut(|a, b| a.original_face_index == b.original_face_index)
        .filter_map(|chunk| match chunk {
            [a, b] => Some([a, b]),
            _ => None,
        })
        .filter(|face| {
            face.iter()
                .all(|v| !v.flags.contains(TriangleFlags::MARK_DEGENERATE))
        })
        .filter(|face| {
            face.iter()
                .filter(|v| v.flags.contains(TriangleFlags::ORIENT_PRESERVING))
                .count()
                == 1
        })
        .for_each(|[mut a, mut b]| {
            // if this happens the quad has extremely bad mapping!!
            if b.flags.contains(TriangleFlags::GROUP_WITH_ANY)
                || calc_tex_area(geometry, &vertices[a.vertex_offset..(a.vertex_offset + 3)])
                    >= calc_tex_area(geometry, &vertices[b.vertex_offset..(b.vertex_offset + 3)])
            {
                (a, b) = (b, a);
            }

            b.flags.set(
                TriangleFlags::ORIENT_PRESERVING,
                a.flags.contains(TriangleFlags::ORIENT_PRESERVING),
            );
        });

    build_neighbors_fast(triangles, vertices);
}

fn build_neighbors_fast<V: VectorSpace>(triangles: &mut [TriangleInfo<V>], vertices: &[usize]) {
    // build array of edges
    let mut edges = vertices
        .iter()
        .enumerate()
        .map(|(a, &i0)| {
            let face = a / 3;
            let b = 3 * face + (a + 1).rem_euclid(3);
            let i1 = vertices[b];

            Edge {
                i0: i0.min(i1),
                i1: i0.max(i1),
                face,
            }
        })
        .collect::<Vec<_>>();

    edges.sort();

    let mut iter = edges.iter().map(|&edge| {
        let vertices = vertices.iter().copied().skip(edge.face * 3).take(3).cycle();

        let (edgenum, (i0, i1)) = vertices
            .clone()
            .zip(vertices.skip(1))
            .enumerate()
            .find(|&(_, (a, b))| (a == edge.i0 || a == edge.i1) && (b == edge.i0 || b == edge.i1))
            .unwrap();

        (edge, edgenum, i0, i1)
    });

    // pair up, adjacent triangles
    while let Some((a, edgenum_a, i0_a, i1_a)) = iter.next() {
        if triangles[a.face].face_neighbors[edgenum_a].is_some() {
            continue;
        }

        // get true index ordering
        let search = iter
            .clone()
            .take_while(|(b, _, _, _)| a.i0 == b.i0 && a.i1 == b.i1)
            .find_map(|(b, edgenum_b, i1_b, i0_b)| {
                let coincident = i0_a == i0_b && i1_a == i1_b;
                let unassigned = triangles[b.face].face_neighbors[edgenum_b].is_none();
                (coincident && unassigned).then_some((edgenum_b, b))
            });

        if let Some((edgenum_b, b)) = search {
            triangles[a.face].face_neighbors[edgenum_a] = Some(b.face);
            triangles[b.face].face_neighbors[edgenum_b] = Some(a.face);
        }
    }
}

/// returns the texture area times 2
fn calc_tex_area<G: Geometry>(geometry: &G, i: &[usize]) -> f32 {
    let t = [i[0], i[1], i[2]].map(|i| geometry.tex_coord_by_index(i).into());
    let a = t[1][0] - t[0][0];
    let b = t[1][1] - t[0][1];
    let c = t[2][0] - t[0][0];
    let d = t[2][1] - t[0][1];
    let signed_area_stx2 = a * d - b * c;
    abs(signed_area_stx2)
}

/// Stores triangle and vertex information split into good and degenerate slices.
struct Info<'a, V: VectorSpace> {
    good_triangles: &'a mut [TriangleInfo<V>],
    degenerate_triangles: &'a mut [TriangleInfo<V>],
    good_vertices: &'a mut [usize],
    degenerate_vertices: &'a mut [usize],
}

/// Sort the provided triangle info and vertex index lists into good and degenerate
/// halves.
/// The returned [`Info`] stores each part as slices of the provided lists.
///
/// Mark all triangle pairs that belong to a quad with only one good triangle.
/// These need special treatment in DegenEpilogue().
/// Additionally, move all good triangles to the start of pTriInfos[] and
/// piTriListIn[] without changing order and put the degenerate triangles last.
fn split_degenerate<'a, G: Geometry>(
    geometry: &G,
    triangles: &'a mut [TriangleInfo<G::Space>],
    vertices: &'a mut [usize],
) -> Info<'a, G::Space> {
    // Mark & count all degenerate triangles
    let degen = triangles
        .iter_mut()
        .zip(vertices.chunks_exact(3))
        .filter(|(_info, i)| {
            let p = [i[0], i[1], i[2]].map(|i| geometry.position_by_index(i));
            p[0] == p[1] || p[0] == p[2] || p[1] == p[2]
        })
        .map(|(info, _)| info.flags.set(TriangleFlags::MARK_DEGENERATE, true))
        .count();

    // locate quads with only one good triangle
    triangles
        .chunk_by_mut(|a, b| a.original_face_index == b.original_face_index)
        .filter_map(|chunk| match chunk {
            [a, b] => Some([a, b]),
            _ => None,
        })
        .filter(|face| {
            face.iter()
                .filter(|v| v.flags.contains(TriangleFlags::MARK_DEGENERATE))
                .count()
                == 1
        })
        .flatten()
        .for_each(|v| v.flags |= TriangleFlags::QUAD_ONE_DEGEN_TRI);

    // reorder list so all degen triangles are moved to the back
    // without reordering the good triangles
    triangles.sort_by(|a, b| {
        match (
            a.flags.contains(TriangleFlags::MARK_DEGENERATE),
            b.flags.contains(TriangleFlags::MARK_DEGENERATE),
        ) {
            (true, true) | (false, false) => core::cmp::Ordering::Equal,
            (true, false) => core::cmp::Ordering::Greater,
            (false, true) => core::cmp::Ordering::Less,
        }
    });

    let (good_triangles, degenerate_triangles) = triangles.split_at_mut(triangles.len() - degen);
    let (good_vertices, degenerate_vertices) = vertices.split_at_mut(3 * good_triangles.len());

    Info {
        good_triangles,
        degenerate_triangles,
        good_vertices,
        degenerate_vertices,
    }
}

/// Make a welded index list of identical positions and attributes `(pos, norm, texc)`.
fn weld_vertices<G: Geometry>(geometry: &G, indices: &mut [usize]) {
    let mut originals = alloc::collections::BTreeMap::new();

    for index in indices {
        let p = geometry
            .position_by_index(*index)
            .into()
            .map(|v| v.to_ne_bytes());
        let n = geometry
            .normal_by_index(*index)
            .into()
            .map(|v| v.to_ne_bytes());
        let t = geometry
            .tex_coord_by_index(*index)
            .into()
            .map(|v| v.to_ne_bytes());

        // texture coordinates are only 2D
        let t = [t[0], t[1]];

        *index = *originals.entry((p, n, t)).or_insert(*index);
    }
}
