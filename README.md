# MonoAlg3D ![build](https://github.com/rsachetto/MonoAlg3D_C/actions/workflows/build.yml/badge.svg)

The MonoAlg3D is a program for solving the 3D monodomain equation by applying the Finite Volume Method.

# Pre-Requisites

  - Linux-based operating system
  - Nvidia Driver 
  - CUDA library

### Setting the enviroment

Ubuntu: Refer to [the ubuntu guide](guide-monoalg3d-ubuntu.md)

Fedora: Refer to [the fedora guide](guide-monoalg3d-fedora.md)

### Compile
```sh
$ ./build.sh
```
The binary files will be saved in the ```bin``` folder.

# Running examples
----
```sh
$ bin/MonoAlg3D -c example_configs/cuboid_ohara.ini 
```

The output will be saved in the VTK format. In order to see the results you can use Paraview (https://www.paraview.org/) or the compiled visualization tool in ```bin/MonoAlg3D_visualizer```. You can also set the output to plain text, by changing the section ```save_result``` in example_configs/cuboid_ohara.ini to:

```ìni
[save_result]
print_rate=250
output_dir=./outputs/tmp_cube
main_function=save_as_text_or_binary
binary=false
```

In the plain text format we have:

- Each line represents a Volume
- Each volume is represented by its center point (X, Y, and Z), the value of half of its side length on x, y and z and the calculated V

Example file:

```
850,850,950,50,50,50, -85
850,950,950,50,50,50, -85
850,950,850,50,50,50, -85
```

This file represents 3 volumes with 100 micrometer of side. The first volume is centered at  at 850,850,950 and the calculated V is -85 mV.

# How to cite:
----

Oliveira RS, Rocha BM, Burgarelli D, Meira Jr W, Constantinides C, dos Santos RW. Performance evaluation of GPU parallelization, space‐time adaptive algorithms, and their combination for simulating cardiac electrophysiology. Int J Numer Meth Biomed Engng. 2018;34:e2913. https://doi.org/10.1002/cnm.2913

# Credits
----
[Heart icons created by phatplus - Flaticon](https://www.flaticon.com/free-icons/heart)
