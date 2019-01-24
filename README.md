# MonoAlg3D
The MonoAlg3D is a program for solving the 3D monodomain equation by applying the Finite Volume Method.

# Pre-Requisites

  - Linux-based operating system
  - Nvidia Driver 
  - CUDA library
  - CMake
  - The Visualization Toolkit library (VTK)
  - Paraview
  - GTK+3 library

### Setting the enviroment

This guide will provide instructions for running the MonoAlg3D code on a Fedora 28 machine.

#### Step 1 - Update your packages

```sh
$ sudo dnf update
```
#### Step 2 - Change graphic mode to GNOME - Xorg

In order to properly install the NVIDIA drivers exit your current session and change the graphic mode to GNOME on Xorg (click on the gear symbol).

#### Step 3 - Install the NVIDIA drivers

##### Add the RPMFusion repositories
```sh
$ sudo dnf install https://download1.rpmfusion.org/free/fedora/rpmfusion-free-release-$(rpm -E %fedora).noarch.rpm https://download1.rpmfusion.org/nonfree/fedora/rpmfusion-nonfree-release-$(rpm -E %fedora).noarch.rpm
```
##### Install the NVIDIA proprietary drivers
```sh
$ sudo dnf install xorg-x11-drv-nvidia akmod-nvidia
```
**Reboot** your machine and check if the installation was successful by going in **Settings > Details**. Your graphics card name must appear as the default Graphic. 

#### Step 4 - Install the CUDA libraries

##### Install the OpenGL libraries
```sh
$ sudo dnf install freeglut-devel
$ sudo dnf install libXt-devel libXmu-devel libXi-devel
```

##### Install the CUDA initial packages
```sh
$ sudo dnf install xorg-x11-drv-nvidia-cuda
```
You can check if everything is going well by running the command ```nvidia-smi```. It should print out some information about your graphic cards. 

##### Download and install the CUDA Toolkit 9.1 for the Fedora 25
```sh
$ wget -nc https://developer.nvidia.com/compute/cuda/9.1/Prod/local_installers/cuda_9.1.85_387.26_linux
```
Go to your Download folder and move the file to a separate folder
```sh
$ mkdir CUDA-build
$ cd Downloads; mv cuda_9.1.85_387.26_linux.run ../CUDA-build 
$ cd CUDA-build
```
Build a temporary folder to store the files from the installation of the Toolkit
```sh
$ sudo mount -o bind /var/tmp /tmp
```
Change the permission of the ```runfile```
```sh
$ chmod +x cuda_9.1.85_387.26_linux.run
```
Install the Toolkit
```sh
$ sudo ./cuda_9.1.85_387.26_linux.run --tmpdir="/tmp" --override
```
When the installer ask you about for the installation of the CUDA driver, select **No**, since you already installed this package.
##### Solve some compatibility problems

###### Edit the ```host_config.h``` file:
Open the file with ```sudo``` and change the line 113 to:
```
#if __GNUC__ > 10 
```
###### Edit the ```/usr/include/bits/floatn.h``` file:
Make a backup of this file.
```sh
$ sudo cp /usr/include/bits/floatn.h /usr/include/bits/floatn-BACKUP.h
```
Then, add the following lines after the line 35:
```sh
#if CUDART_VERSION
#undef __HAVE_FLOAT128
#define __HAVE_FLOAT128 0
#endif
```
###### Edit your ```.bash_profile``` file:
Add the following lines:
```sh
PATH=$PATH:$HOME/.local/bin:$HOME/bin:/usr/local/cuda-9.1/bin
LIBRARY_PATH=$LIBRARY_PATH:$HOME/lib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-9.1/lib64:$HOME/lib
CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$HOME/NVIDIA_CUDA-9.1_Samples/common/inc:$HOME/include
 
export PATH
export LIBRARY_PATH
export LD_LIBRARY_PATH
export CPLUS_INCLUDE_PATH
```
###### Install the GCC 5.3.0:
To avoid some incompatibilities with the new version of GCC we need to use a previous version of GNU/GCC compiler. Go to this URL and download and install the RPM package:

**GNU/GCC compiler v5.3.0**
https://drive.google.com/file/d/0B7S255p3kFXNbTBneHgwSzBodFE/view

###### Make a symbolic link for the CUDA compiler:
```sh
$ sudo ln -s /usr/bin/gcc53 /usr/local/cuda/bin/gcc
$ sudo ln -s /usr/bin/g++53 /usr/local/cuda/bin/g++
```

**Reboot** your machine.

#### Step 5 - Install CMake
To install the CMake package just run the following command:
```sh
$ sudo dnf install cmake
```

#### Step 6 - Install VTK libraries

##### Download the version 7.1.1 of VTK:
```sh
$ wget -nc https://www.vtk.org/files/release/7.1/VTK-7.1.1.tar.gz
```
##### Extract the file in your ```home``` directory:
```sh
$ cd; mv Downloads/VTK-7.1.1.tar.gz ~; tar -xf VTK-7.1.1.tar.gz  
```
##### Make a new folder to store the build files:
```sh
$ mkdir VTK-build; cd VTK-build  
```
##### Build the VTK files:
```sh
$ ccmake ../VTK-7.1.1 
```
Press the 'c' key to configure. Next apply the configuration settings by hitting 'c' again. Generate the Makefile by pressing 'g'.

Finally, we can build the VTK files by running:
```sh
$ sudo make install 
```
##### Optional
You can speedup the building process by adding the ```-j``` flag to the ```make``` command. For example, running the below command will build the files using 3 cores from your CPU in parallel.
```sh
$ sudo make -j3 install 
```
#### Step 7 - Install Paraview
The Paraview program will allow the visualization of the simulations, to installed simply execute:
```sh
$ sudo dnf install paraview
```
#### Step 8 - Install GTK+3 library
For the GUI interface we need these libraries.
```sh
$ sudo dnf install gtk3-devel gtksourceview3-devel
```

# Compile
----
Change the current version of GCC to 5.3.0
```sh
$ export CC=/usr/bin/gcc53
$ export CXX=/usr/bin/g++53
```

Build the binary
```sh
$ mkdir build
$ cd build
$ cmake ..
$ make
```
The binary file for execution will be saved in the ```bin``` folder.

# Running examples
----
```sh
$ bin/MonoAlg3D -c example_configs/cuboid_ohara.ini 
```

The output will be saved in the VTK format. In order to see the results you can use Paraview (https://www.paraview.org/). You can also set the output to plain text, by changing the option ```vtk_output``` to false in the configuration file. The text format is defined as following:

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
----

Oliveira RS, Rocha BM, Burgarelli D, Meira Jr W, Constantinides C, dos Santos RW. Performance evaluation of GPU parallelization, space‐time adaptive algorithms, and their combination for simulating cardiac electrophysiology. Int J Numer Meth Biomed Engng. 2018;34:e2913. https://doi.org/10.1002/cnm.2913


