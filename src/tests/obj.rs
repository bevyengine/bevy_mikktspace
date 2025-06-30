use glam::Vec3;
use wavefront_obj::obj::{ObjSet, Primitive, parse};

struct GlamSpace;

impl crate::VectorSpace for GlamSpace {
    type Vec3 = Vec3;

    fn length(a: Self::Vec3) -> f32 {
        a.length()
    }

    fn angle_between(a: Self::Vec3, b: Self::Vec3) -> f32 {
        let a = Self::normalize_or_zero(a);
        let b = Self::normalize_or_zero(b);
        let cos = Self::dot(a, b).clamp(-1.0, 1.0);
        (cos as f64).acos() as f32
    }
}

struct WavefrontObject(wavefront_obj::obj::Object);

impl crate::Geometry for WavefrontObject {
    type Space = GlamSpace;

    fn num_faces(&self) -> usize {
        self.0.geometry.iter().flat_map(|g| g.shapes.iter()).count()
    }

    fn num_vertices_of_face(&self, face: usize) -> usize {
        let face = self
            .0
            .geometry
            .iter()
            .flat_map(|g| g.shapes.iter())
            .nth(face)
            .unwrap();
        match face.primitive {
            Primitive::Point(_) => 1,
            Primitive::Line(_, _) => 2,
            Primitive::Triangle(_, _, _) => 3,
        }
    }

    fn position(&self, face: usize, vert: usize) -> Vec3 {
        let shape = self
            .0
            .geometry
            .iter()
            .flat_map(|g| g.shapes.iter())
            .nth(face)
            .unwrap();

        let indices = match shape.primitive {
            Primitive::Point(a) => &[a][..],
            Primitive::Line(a, b) => &[a, b][..],
            Primitive::Triangle(a, b, c) => &[a, b, c][..],
        };

        let index = indices.iter().map(|(v, _t, _n)| *v).nth(vert).unwrap();

        let vertex = self.0.vertices[index];

        [vertex.x as f32, vertex.y as f32, vertex.z as f32].into()
    }

    fn normal(&self, face: usize, vert: usize) -> Vec3 {
        let shape = self
            .0
            .geometry
            .iter()
            .flat_map(|g| g.shapes.iter())
            .nth(face)
            .unwrap();

        let indices = match shape.primitive {
            Primitive::Point(a) => &[a][..],
            Primitive::Line(a, b) => &[a, b][..],
            Primitive::Triangle(a, b, c) => &[a, b, c][..],
        };

        let index = indices
            .iter()
            .map(|(_v, _t, n)| *n)
            .nth(vert)
            .unwrap()
            .unwrap();

        let normal = self.0.normals[index];

        [normal.x as f32, normal.y as f32, normal.z as f32].into()
    }

    fn tex_coord(&self, face: usize, vert: usize) -> [f32; 2] {
        let shape = self
            .0
            .geometry
            .iter()
            .flat_map(|g| g.shapes.iter())
            .nth(face)
            .unwrap();

        let indices = match shape.primitive {
            Primitive::Point(a) => &[a][..],
            Primitive::Line(a, b) => &[a, b][..],
            Primitive::Triangle(a, b, c) => &[a, b, c][..],
        };

        let index = indices
            .iter()
            .map(|(_v, t, _n)| *t)
            .nth(vert)
            .unwrap()
            .unwrap();

        let tvertex = self.0.tex_vertices[index];

        [tvertex.u as f32, tvertex.v as f32]
    }
}

#[test]
fn cube() {
    let file = include_str!("../../data/cube.obj");
    let ObjSet { objects, .. } = parse(file).expect("must be able to parse sample data");
    for object in objects {
        let mut object = WavefrontObject(object);
        crate::generate_tangents(&mut object);
    }
}
