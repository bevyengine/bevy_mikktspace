#![no_main]

use libfuzzer_sys::fuzz_target;
use bevy_mikktspace_fuzz::Geometry;

fuzz_target!(|input: (Geometry, f32)| {
    let (value, angular_threshold) = input;
    let angular_threshold = angular_threshold % 180.;

    let reference = {
        let mut value = value.clone();
        mikktspace_sys::gen_tang_space(&mut value, angular_threshold);
        value
    };

    let value = {
        let mut value = value;
        let linear_threshold = angular_threshold.to_radians().cos();
        let _ = bevy_mikktspace::generate_tangents_with_threshold(&mut value, linear_threshold);
        value
    };

    reference.assert_eq(&value);
});
