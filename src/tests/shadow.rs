//! Provides [`Shadow`], which wraps a [`Geometry`] implementation and will compare
//! results to [`mikktspace_sys`].

use alloc::collections::BTreeMap;

use crate::Geometry;

/// Wrapper that implements [`mikktspace_sys::MikkTSpaceInterface`] and [`Geometry`],
/// but will internally store results from [`mikktspace_sys`] calculations to compare
/// against this crate's implementation.
pub(crate) struct Shadow<'a, G: Geometry> {
    inner: &'a mut G,
    /// Results computed using [`mikktspace_sys`] to use as a reference.
    reference_results: BTreeMap<Key, Entry>,
}

impl<'a, G: Geometry> Shadow<'a, G> {
    /// Construct a new [`Shadow`] and immediately compute a reference set of results.
    pub(crate) fn new(inner: &'a mut G) -> Self {
        let mut this = Self {
            inner,
            reference_results: BTreeMap::new(),
        };
        mikktspace_sys::gen_tang_space_default(&mut this);
        this
    }
}

#[derive(Clone, Copy, PartialEq, Debug)]
struct Entry {
    tangent: [f32; 3],
    bi_tangent: [f32; 3],
    f_mag_s: f32,
    f_mag_t: f32,
    bi_tangent_preserves_orientation: bool,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug, PartialOrd, Ord)]
struct Key {
    face: usize,
    vert: usize,
}

impl<G: Geometry> mikktspace_sys::MikkTSpaceInterface for Shadow<'_, G> {
    fn get_num_faces(&self) -> usize {
        self.inner.num_faces()
    }

    fn get_num_vertices_of_face(&self, face: usize) -> usize {
        self.inner.num_vertices_of_face(face)
    }

    fn get_position(&self, face: usize, vert: usize) -> [f32; 3] {
        self.inner.position(face, vert).into()
    }

    fn get_normal(&self, face: usize, vert: usize) -> [f32; 3] {
        self.inner.normal(face, vert).into()
    }

    fn get_tex_coord(&self, face: usize, vert: usize) -> [f32; 2] {
        self.inner.tex_coord(face, vert)
    }

    fn set_tspace(
        &mut self,
        tangent: [f32; 3],
        bi_tangent: [f32; 3],
        mag_s: f32,
        mag_t: f32,
        is_orientation_preserving: bool,
        face: usize,
        vert: usize,
    ) {
        let key = Key { face, vert };
        let entry = Entry {
            tangent,
            bi_tangent,
            f_mag_s: mag_s,
            f_mag_t: mag_t,
            bi_tangent_preserves_orientation: is_orientation_preserving,
        };
        self.reference_results.insert(key, entry);
        self.inner.set_tangent(
            tangent.into(),
            bi_tangent.into(),
            mag_s,
            mag_t,
            is_orientation_preserving,
            face,
            vert,
        );
    }
}

impl<G: Geometry> Geometry for Shadow<'_, G> {
    type Space = G::Space;

    fn num_faces(&self) -> usize {
        self.inner.num_faces()
    }

    fn num_vertices_of_face(&self, face: usize) -> usize {
        self.inner.num_vertices_of_face(face)
    }

    fn position(&self, face: usize, vert: usize) -> <G::Space as crate::VectorSpace>::Vec3 {
        self.inner.position(face, vert)
    }

    fn normal(&self, face: usize, vert: usize) -> <G::Space as crate::VectorSpace>::Vec3 {
        self.inner.normal(face, vert)
    }

    fn tex_coord(&self, face: usize, vert: usize) -> [f32; 2] {
        self.inner.tex_coord(face, vert)
    }

    fn set_tangent(
        &mut self,
        tangent: <G::Space as crate::VectorSpace>::Vec3,
        bi_tangent: <G::Space as crate::VectorSpace>::Vec3,
        f_mag_s: f32,
        f_mag_t: f32,
        bi_tangent_preserves_orientation: bool,
        face: usize,
        vert: usize,
    ) {
        let key = Key { face, vert };
        let entry = Entry {
            tangent: tangent.into(),
            bi_tangent: bi_tangent.into(),
            f_mag_s,
            f_mag_t,
            bi_tangent_preserves_orientation,
        };

        let Some(reference) = self.reference_results.get(&key) else {
            panic!("A reference tangent was not calculated for {key:?}");
        };

        assert_eq!(reference, &entry);

        self.inner.set_tangent(
            tangent,
            bi_tangent,
            f_mag_s,
            f_mag_t,
            bi_tangent_preserves_orientation,
            face,
            vert,
        );
    }
}
