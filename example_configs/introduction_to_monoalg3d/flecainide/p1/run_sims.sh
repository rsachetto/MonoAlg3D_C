#!/bin/bash
files=(
    example_configs/introduction_to_monoalg3d/flecainide/p1/EX00_locally_applied_flecainide_0_1_p1.ini
    example_configs/introduction_to_monoalg3d/flecainide/p1/EX00_locally_applied_flecainide_1_p1.ini
    example_configs/introduction_to_monoalg3d/flecainide/p1/EX00_locally_applied_flecainide_75_p1.ini
    example_configs/introduction_to_monoalg3d/flecainide/p1/EX00_locally_applied_flecainide_reference.ini
)

for file in "${files[@]}"; do
    if [[ -f $file ]]; then
        ./bin/MonoAlg3D -c "$file"
    fi
done
