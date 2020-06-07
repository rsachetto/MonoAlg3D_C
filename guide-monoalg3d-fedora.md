### GUIDE TO CONFIGURE THE MONOALG_3D ON A FEDORA 28 MACHINE

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

----
Change the current version of GCC to 5.3 by changing the following variables in bsbash/build_functions.sh (For manjaro no changes are needed. For Fedora you need the following changes.)
```sh
$ C_COMPILER=/usr/bin/gcc53
$ CXX_COMPILER=/usr/bin/g++53
```