#![no_main]

use libfuzzer_sys::fuzz_target;
use bevy_mikktspace_fuzz::Geometry;

fuzz_target!(|value: Geometry| {
    let reference = {
        let mut value = value.clone();
        mikktspace_sys::gen_tang_space_default(&mut value);
        value
    };

    let value = {
        let mut value = value;
        let _ = bevy_mikktspace::generate_tangents(&mut value);
        value
    };

    reference.assert_eq(&value);
});
