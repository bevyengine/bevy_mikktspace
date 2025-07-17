use arbitrary::{Arbitrary, Unstructured};

#[derive(Debug, Clone, Copy, Arbitrary)]
pub struct Vertex {
    pub position: [f32; 3],
    pub normal: [f32; 3],
    pub tex_coord: [f32; 2],
    pub tangent: [f32; 3],
}

impl PartialEq for Vertex {
    fn eq(&self, other: &Vertex) -> bool {
        if self.position.map(|p| p.to_ne_bytes()) != other.position.map(|p| p.to_ne_bytes()) {
            return false;
        }

        if self.normal.map(|p| p.to_ne_bytes()) != other.normal.map(|p| p.to_ne_bytes()) {
            return false;
        }

        if self.tex_coord.map(|p| p.to_ne_bytes()) != other.tex_coord.map(|p| p.to_ne_bytes()) {
            return false;
        }

        self.tangent.map(|p| p.to_ne_bytes()) == other.tangent.map(|p| p.to_ne_bytes())
    }
}

#[derive(Debug, Clone, PartialEq, Arbitrary)]
pub enum Face {
    Triangle([usize; 3]),
    Quad([usize; 4]),
    Arbitrary(Vec<usize>),
}

impl core::ops::Deref for Face {
    type Target = [usize];

    fn deref(&self) -> &Self::Target {
        match self {
            Self::Triangle(inner) => &*inner,
            Self::Quad(inner) => &*inner,
            Self::Arbitrary(inner) => inner.as_slice(),
        }
    }
}

impl core::ops::DerefMut for Face {
    fn deref_mut(&mut self) -> &mut Self::Target {
        match self {
            Self::Triangle(inner) => &mut *inner,
            Self::Quad(inner) => &mut *inner,
            Self::Arbitrary(inner) => inner.deref_mut(),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Geometry {
    pub vertices: Vec<Vertex>,
    pub faces: Vec<Face>,
}

impl Geometry {
    pub fn validate(&mut self) -> Result<(), arbitrary::Error> {
        let Self {
            vertices,
            faces,
        } = self;

        // Known failure: no vertices
        if vertices.is_empty() {
            return Err(arbitrary::Error::IncorrectFormat);
        }

        // Known failure: no faces
        if faces.is_empty() {
            return Err(arbitrary::Error::IncorrectFormat);
        }

        // Known failure: face vertex indices out of range
        for face in faces.iter_mut() {
            if face.is_empty() {
                return Err(arbitrary::Error::IncorrectFormat);
            }
            
            for vertex in face.iter_mut() {
                *vertex %= vertices.len();
            }
        }

        // Known failure: NaN values
        if vertices
            .iter()
            .flat_map(|vertex| {
                vertex
                    .position
                    .iter()
                    .copied()
                    .chain(vertex.normal)
                    .chain(vertex.tex_coord)
            })
            .any(|v| v.is_nan())
        {
            return Err(arbitrary::Error::IncorrectFormat);
        }

        // Known failure: Identical positions within a face
        for face in faces.iter() {
            let mut iter = face.iter();
            while let Some(a) = iter.next() {
                for b in iter.clone() {
                    if vertices[*a].position == vertices[*b].position {
                        return Err(arbitrary::Error::IncorrectFormat);
                    }
                }
            }
        }

        // Known failure: non-normal values
        if vertices
            .iter()
            .flat_map(|vertex| {
                vertex
                    .position
                    .iter()
                    .copied()
                    .chain(vertex.normal)
                    .chain(vertex.tex_coord)
            })
            .any(|v| !v.is_normal())
        {
            return Err(arbitrary::Error::IncorrectFormat);
        }

        Ok(())
    }

    /// Compares two geometries totally before panicking _if_ they differ.
    pub fn assert_eq(&self, other: &Self) {
        use core::fmt::Write;

        let mut errors = String::new();

        if self.faces.len() != other.faces.len() {
            let _ = writeln!(&mut errors, "- Expected {} triangles; found {}", self.faces.len(), other.faces.len());
        }

        for (a, b) in self.faces.iter().zip(other.faces.iter()).filter(|(a, b)| a != b) {
            let _ = writeln!(&mut errors, "- Expected {:?} triangle; found {:?}", a, b);
        }

        if self.vertices.len() != other.vertices.len() {
            let _ = writeln!(&mut errors, "- Expected {} vertices; found {}", self.vertices.len(), other.vertices.len());
        }

        for (a, b) in self.vertices.iter().zip(other.vertices.iter()).filter(|(a, b)| a != b) {
            let _ = writeln!(&mut errors, "  - Differing vertex:");
            if a.position != b.position {
                let _ = writeln!(&mut errors, "    - Expected position {:?}; found {:?}", a.position, b.position);
            }
            if a.normal != b.normal {
                let _ = writeln!(&mut errors, "    - Expected normal {:?}; found {:?}", a.normal, b.normal);
            }
            if a.tex_coord != b.tex_coord {
                let _ = writeln!(&mut errors, "    - Expected texture coordinate {:?}; found {:?}", a.tex_coord, b.tex_coord);
            }
            if a.tangent != b.tangent {
                let _ = writeln!(&mut errors, "    - Expected tangent {:?}; found {:?}", a.tangent, b.tangent);
            }
        }

        if !errors.is_empty() {
            panic!("Difference Summary:\n{}", errors);
        }
    }
}

impl Arbitrary<'_> for Geometry {
    fn arbitrary(u: &mut Unstructured<'_>) -> Result<Self, arbitrary::Error> {
        let mut value = Self {
            vertices: Vec::<Vertex>::arbitrary(u)?,
            faces: Vec::<Face>::arbitrary(u)?,
        };

        value.validate()?;

        Ok(value)
    }
}

impl mikktspace_sys::MikkTSpaceInterface for Geometry {
    fn get_num_faces(&self) -> usize {
        self.faces.len()
    }

    fn get_num_vertices_of_face(&self, face: usize) -> usize {
        self.faces[face].len()
    }

    fn get_position(&self, face: usize, vert: usize) -> [f32; 3] {
        self.vertices[self.faces[face][vert]].position
    }

    fn get_normal(&self, face: usize, vert: usize) -> [f32; 3] {
        self.vertices[self.faces[face][vert]].normal
    }

    fn get_tex_coord(&self, face: usize, vert: usize) -> [f32; 2] {
        self.vertices[self.faces[face][vert]].tex_coord
    }

    fn set_tspace(
        &mut self,
        tangent: [f32; 3],
        _bi_tangent: [f32; 3],
        _mag_s: f32,
        _mag_t: f32,
        _is_orientation_preserving: bool,
        face: usize,
        vert: usize,
    ) {
        self.vertices[self.faces[face][vert]].tangent = tangent;
    }
}

impl bevy_mikktspace::Geometry for Geometry {
    fn num_faces(&self) -> usize {
        self.faces.len()
    }

    fn num_vertices_of_face(&self, face: usize) -> usize {
        self.faces[face].len()
    }

    fn position(&self, face: usize, vert: usize) -> [f32; 3] {
        self.vertices[self.faces[face][vert]].position
    }

    fn normal(&self, face: usize, vert: usize) -> [f32; 3] {
        self.vertices[self.faces[face][vert]].normal
    }

    fn tex_coord(&self, face: usize, vert: usize) -> [f32; 2] {
        self.vertices[self.faces[face][vert]].tex_coord
    }

    fn set_tangent(
        &mut self,
        tangent_space: Option<bevy_mikktspace::TangentSpace>,
        face: usize,
        vert: usize,
    ) {
        let tangent_space = tangent_space.unwrap_or_default();
        self.vertices[self.faces[face][vert]].tangent = tangent_space.tangent();
    }
}
