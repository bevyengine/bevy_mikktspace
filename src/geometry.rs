use crate::math::VectorSpace;

/// Interface for the mikktspace algorithm to query information about your geometry.
pub trait Geometry {
    type Space: VectorSpace;

    /// Returns the number of faces.
    fn num_faces(&self) -> usize;

    /// Returns the number of vertices of a face.
    fn num_vertices_of_face(&self, face: usize) -> usize;

    /// Returns the position of a vertex.
    fn position(&self, face: usize, vert: usize) -> <Self::Space as VectorSpace>::Vec3;

    /// Returns the normal of a vertex.
    fn normal(&self, face: usize, vert: usize) -> <Self::Space as VectorSpace>::Vec3;

    /// Returns the texture coordinate of a vertex.
    fn tex_coord(&self, face: usize, vert: usize) -> [f32; 2];

    /// Sets the generated tangent for a vertex.
    /// Leave this function unimplemented if you are implementing
    /// `set_tangent_encoded`.
    #[expect(clippy::too_many_arguments, reason = "a lot of arguments are needed")]
    fn set_tangent(
        &mut self,
        tangent: <Self::Space as VectorSpace>::Vec3,
        _bi_tangent: <Self::Space as VectorSpace>::Vec3,
        _f_mag_s: f32,
        _f_mag_t: f32,
        bi_tangent_preserves_orientation: bool,
        face: usize,
        vert: usize,
    ) {
        let tangent = tangent.into();
        let sign = if bi_tangent_preserves_orientation {
            1.0
        } else {
            -1.0
        };
        self.set_tangent_encoded([tangent[0], tangent[1], tangent[2], sign], face, vert);
    }

    /// Sets the generated tangent for a vertex with its bi-tangent encoded as the 'W' (4th)
    /// component in the tangent. The 'W' component marks if the bi-tangent is flipped. This
    /// is called by the default implementation of `set_tangent`; therefore, this function will
    /// not be called by the crate unless `set_tangent` is unimplemented.
    fn set_tangent_encoded(&mut self, _tangent: [f32; 4], _face: usize, _vert: usize) {}
}

pub(crate) trait GeometryExt: Geometry {
    /// Returns the position of a vertex on a face using the provided index.
    fn position_by_index(&self, index: usize) -> <Self::Space as VectorSpace>::Vec3 {
        let (face, vert) = index_to_face_vertex(index);
        self.position(face, vert)
    }

    /// Returns the texture coordinate of a vertex on a face using the provided index.
    fn tex_coord_by_index(&self, index: usize) -> <Self::Space as VectorSpace>::Vec3 {
        let (face, vert) = index_to_face_vertex(index);
        let tex_coord = self.tex_coord(face, vert);
        [tex_coord[0], tex_coord[1], 1.0].into()
    }

    /// Returns the normal of a vertex on a face using the provided index.
    fn normal_by_index(&self, index: usize) -> <Self::Space as VectorSpace>::Vec3 {
        let (face, vert) = index_to_face_vertex(index);
        self.normal(face, vert)
    }

    /// Iterate over all faces by index.
    fn face_indices(&self) -> impl ExactSizeIterator<Item = usize> + Clone {
        0..self.num_faces()
    }

    /// Iterate over all vertices on a provided face.
    fn vertex_indices(&self, face: usize) -> impl ExactSizeIterator<Item = usize> + Clone {
        (0..self.num_vertices_of_face(face)).map(move |i| face_vertex_to_index(face, i))
    }
}

impl<G: Geometry> GeometryExt for G {}

// Mikktspace uses indices internally to refer to and identify vertices, these utility functions
// make it easier to work with these indices.

/// Generate a vertex index for the Nth vertex of the Nth face.
fn face_vertex_to_index(face: usize, vertex: usize) -> usize {
    face << 2 | vertex & 0x3
}

/// Reverse of `face_vertex_to_index`.
fn index_to_face_vertex(index: usize) -> (usize, usize) {
    (index >> 2, index & 0x3)
}
