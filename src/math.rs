use core::ops::{Add, Mul, Neg, Sub};

pub trait VectorSpace {
    type Vec3: Clone
        + Copy
        + Neg<Output = Self::Vec3>
        + Add<Self::Vec3, Output = Self::Vec3>
        + Sub<Self::Vec3, Output = Self::Vec3>
        + Mul<Self::Vec3, Output = Self::Vec3>
        + PartialEq
        + From<[f32; 3]>
        + Into<[f32; 3]>;

    fn length(a: Self::Vec3) -> f32;

    fn angle_between(a: Self::Vec3, b: Self::Vec3) -> f32;

    fn scale(this: Self::Vec3, t: f32) -> Self::Vec3 {
        this.into().map(|x| x * t).into()
    }

    fn dot(this: Self::Vec3, rhs: Self::Vec3) -> f32 {
        let [ax, ay, az] = this.into();
        let [bx, by, bz] = rhs.into();
        ax * bx + ay * by + az * bz
    }

    fn length_squared(a: Self::Vec3) -> f32 {
        Self::dot(a, a)
    }

    fn distance_squared(a: Self::Vec3, b: Self::Vec3) -> f32 {
        Self::length_squared(a - b)
    }

    fn project(a: Self::Vec3, b: Self::Vec3) -> Self::Vec3 {
        a - Self::scale(b, Self::dot(b, a))
    }

    fn normalize_or_zero(a: Self::Vec3) -> Self::Vec3 {
        let [x, y, z] = a.into();
        if not_zero(x) || not_zero(y) || not_zero(z) {
            Self::scale(a, Self::length(a).recip())
        } else {
            a
        }
    }
}

pub(crate) fn abs(value: f32) -> f32 {
    if value < 0.0 { -value } else { value }
}

pub(crate) fn not_zero(fx: f32) -> bool {
    abs(fx) > 1.1754944e-38f32
}
