# Compile

$ mkdir build

$ cd build

$ cmake ../

$ make

The binary program will be generated at bin folder.

# Running example:

bin/MonoAlg3D -c example_configs/cuboid_ohara.ini

The output will be saved in the VTK format. In order to see the results you can use Paraview (https://www.paraview.org/). You can also set the output to plain text, by changing the option vtk_output to false in the configuration file. The text format is defined as follwing:

- Each line represents a Volume
- Each volume is represented by its center point (X, Y, and Z), the value of halt of its side length and the calculated V

Example file:

```
850,850,950,50,-85
850,950,950,50,-85
850,950,850,50,-85
```

This file represents 3 volumes with 100 micrometer of side. The first volume is centered at  at 850,850,950 and the calculated V is -85 mV.

# How to cite:

Oliveira RS, Rocha BM, Burgarelli D, Meira Jr W, Constantinides C, dos Santos RW. Performance evaluation of GPU parallelization, space‐time adaptive algorithms, and their combination for simulating cardiac electrophysiology. Int J Numer Meth Biomed Engng. 2018;34:e2913. https://doi.org/10.1002/cnm.2913
