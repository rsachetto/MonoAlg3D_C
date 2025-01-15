## GUIDE TO CONFIGURE MONOALG_3D ON A FEDORA 40 MACHINE

----------------------------------------------------------

This guide provides instructions to configure the MonoAlg3D solver in a fresh instalation of Fedora 40 with all its features working.

## 1) Download and build a USB boot drive

### Download the Fedora 40 ISO image for your system:

- https://dl.fedoraproject.org/pub/fedora/linux/releases/40/Workstation/

- Burn the ISO image in a USB device using the Rufus program (Windows) or any other similar program.

- Install Fedora 40 on the target machine

### Perform a system update:

```sh
$ sudo dnf update
```

### Change the graphic mode to "GNOME over Xorg"

- Close your current session
- Click on your user
- Click on the gear symbol and select the _"GNOME over Xorg"_
- Begin your session again

### Reboot your system

```sh
$ reboot
```

### Install the GNOME-Tweak-Tools

```sh
$ sudo dnf install gnome-tweak-tool
```

### Open Configurations

- Disable energy saving options
- Disable the suspend screen when notebook is closed

### Install gedit
```sh
$ sudo dnf gedit
```

### Install the Resources program to replace the System Monitor

- Search for the 'Programs' software
- Open 'Programs' and search for the Resources program
- Install using the Flathub package

### Install the Engrampa package manager

- Get the 'Engrampa' software from Programs

## 2) Install RPMFusion

### Add the packages to your system:

```sh
$ sudo dnf install https://download1.rpmfusion.org/free/fedora/rpmfusion-free-release-$(rpm -E %fedora).noarch.rpm https://download1.rpmfusion.org/nonfree/fedora/rpmfusion-nonfree-release-$(rpm -E %fedora).noarch.rpm
$ sudo dnf update
```

- Reboot your system

```sh
$ reboot
```

- Enter your BIOS by pressing F12 (check your computer key)
- Disable the Secure Boot
- Save your changes and exit

## 3) Install the NVIDIA driver

### Run the following commands

```sh
$ sudo dnf update
$ sudo dnf install akmod-nvidia
$ sudo dnf install xorg-x11-drv-nvidia-cuda
```

### Wait around 5 to 8 minutes in order to the NVIDIA kernels get built

- Check if the kernels are ready

```sh
$ modinfo -F version nvidia
```

- This should return the version of the driver such as 465.27

### Reboot your system

```sh
$ reboot
```

### Check if the NVIDIA driver is correctly installed 

- Go to Settings > System > About > System Details
- Your NVIDIA Graphics Card should appear in Graphics field.

## 4) Install CUDA libraries

### Download the CUDA runfile for your operating system

```sh
$ cd ~; mkdir CUDA-Install; cd CUDA-Install
$ wget -nc https://developer.download.nvidia.com/compute/cuda/12.6.3/local_installers/cuda_12.6.3_560.35.05_linux.run
```

- Here, I am using the CUDA version 12.6.3.

### Install dependencies

```sh
$ sudo dnf install gcc-c++ mesa-libGLU-devel libX11-devel libXi-devel libXmu-devel 
$ sudo dnf install freeglut freeglut-devel
$ sudo dnf install gcc13-c++-13.3.1-2.fc40.x86_64
```

### Run the CUDA runfile

```sh
$ chmod +x cuda_12.6.3_560.35.05_linux.run
$ sudo ./cuda_12.6.3_560.35.05_linux.run
```

### Configure the installation

- Accept the terms
- **IMPORTANT!** Disable the NVIDIA Driver installation
- Enable the installation: CUDA Toolkit, CUDA Samples, CUDA Demo, CUDA Documentation
- Install the package

### Post installation steps

#### Log as root user

```sh
$ su -
```

- If is the first time logging as root, you can set a password to the _"root"_ user with:

```sh
$ sudo passwd root
```

#### Set the PATH and LD_LIBRARY_PATH by creating a shell configuration file for CUDA

```sh
cat << EOF > /etc/profile.d/cuda.sh
pathmunge /usr/local/cuda/bin before

if [ -z "${LD_LIBRARY_PATH}" ]; then
    LD_LIBRARY_PATH=/usr/local/cuda/lib64
else
    LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH
fi

export PATH LD_LIBRARY_PATH
EOF
```

#### Logout and login again. Check the variables

```sh
$ echo $PATH
$ echo $LD_LIBRARY_PATH
```

- The new paths must appear.

#### Test the CUDA library

```sh
$ cd /home/<username>/NVIDIA_CUDA-12.6_Samples/1_Utilities/deviceQuery
$ make
$ ./deviceQuery
```

- This should print some information about your GPU hardware.

## 5) Configure MonoAlg3D_C

### Download remaining dependencies

```sh
$ sudo dnf install libXcursor-devel libXrandr-devel libXinerama-devel
```

### Make a symbolic link for GCC version 13

```sh
$ sudo ln -s /usr/bin/gcc-13 /usr/local/cuda/bin/gcc
$ sudo ln -s /usr/bin/g++-13 /usr/local/cuda/bin/g++
```

### Build and compile MonoAlg3D

```sh
$ ./build.sh -f
```

**IMPORTANT!** If you stop here you already have a working MonoAlg3D binary executable in the _/bin_ folder

### [OPTIONAL] Install and configure MPI for batch simulations

#### Download OpenMPI

```sh
$ cd ~; mkdir OpenMPI-build; cd OpenMPI-build
$ wget -nc https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.6.tar.gz
```

#### Build OpenMPI

```sh
$ gunzip -c openmpi-5.0.6.tar.gz | tar xf -
$ cd openmpi-5.0.6
$ ./configure --prefix=/usr/local
$ sudo make all install
```

#### Configure OpenMPI for MonoAlg3D

- Change the following function in the _'bsbash/find_functions.sh'_ file:

```sh
FIND_MPI () {
	MPI_FOUND="y"
	MPI_LIBRARIES="mpi"
	MPI_LIBRARY_PATH=$OMPI_LIBRARY_PATH
	MPI_INCLUDE_PATH=$OMPI_INCLUDE_PATH
}
```

- Rebuild MonoAlg3D

```sh
./build.sh -f
```

- Check if you have the _"MonoAlg3D_batch"_ executable inside the _"/bin"_ folder.