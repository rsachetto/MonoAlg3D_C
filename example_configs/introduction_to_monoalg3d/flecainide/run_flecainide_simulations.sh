#!/bin/bash

# run from MonoAlg3d base directory
# This script runs the simulations for the flecainide example

files=(
  ./example_configs/introduction_to_monoalg3d/flecainide/p1/EX00_locally_applied_flecainide_0_1_p1.ini
  ./example_configs/introduction_to_monoalg3d/flecainide/p1/EX00_locally_applied_flecainide_1_p1.ini
  ./example_configs/introduction_to_monoalg3d/flecainide/p1/EX00_locally_applied_flecainide_75_p1.ini
  #./example_configs/introduction_to_monoalg3d/flecainide/p2/EX00_locally_applied_flecainide_0_1_p2.ini
  #./example_configs/introduction_to_monoalg3d/flecainide/p2/EX00_locally_applied_flecainide_1_p2.ini
  #./example_configs/introduction_to_monoalg3d/flecainide/p2/EX00_locally_applied_flecainide_75_p2.ini
  #./example_configs/introduction_to_monoalg3d/flecainide/p3/EX00_locally_applied_flecainide_0_1_p3.ini
  #./example_configs/introduction_to_monoalg3d/flecainide/p3/EX00_locally_applied_flecainide_1_p3.ini
  #./example_configs/introduction_to_monoalg3d/flecainide/p3/EX00_locally_applied_flecainide_75_p3.ini
)

for file in "${files[@]}"
do
  ./bin/MonoAlg3D -c $file
done