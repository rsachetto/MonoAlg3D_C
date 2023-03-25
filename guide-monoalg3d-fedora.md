## GUIDE TO CONFIGURE THE MONOALG_3D ON A FEDORA 33 MACHINE

----------------------------------------------------------

This guide provides instructions to configure the MonoAlg3D solver in a fresh instalation of Fedora 33 with all its features working.

## 1) Download and build a USB boot drive

### Download the Fedora 33 ISO image:

- https://getfedora.org/pt_BR/

### Burn the ISO image in a USB device using the Rufus program (Windows) or any other program.

### Install Fedora 33 in the target machine

### Perform a system update:

```sh
$ sudo dnf update
```

### Change the graphic mode to "GNOME over Xorg"

- Close your current session
- Click on your user
- Click on the gear symbol and select the "GNOME over Xorg"
- Begin your session again

### Reboot your system

```sh
$ reboot
```

### Install the GNOME-Tweak-Tools

```sh
$ sudo dnf install gnome-tweak-tool
```

### Open the Gnome-Tweak-Tool

- Disable animations to save memory
- Disable the suspend screen when notebook is closed

## 2) Install RPMFusion

### Add the packages to your system:

```sh
$ sudo dnf install https://download1.rpmfusion.org/free/fedora/rpmfusion-free-release-$(rpm -E %fedora).noarch.rpm https://download1.rpmfusion.org/nonfree/fedora/rpmfusion-nonfree-release-$(rpm -E %fedora).noarch.rpm
$ sudo dnf update
```

### Open the Gnome-Programs

- Go to Program Repositories
- Enable the RPMFusion repositories for NVIDIA
- Update the repositories again

```sh
$ sudo dnf update
```

- Reboot your system

```sh
$ reboot
```

## 3) Install the NVIDIA driver

### Reboot your system

- Enter your BIOS by pressing F12 (check your computer key)
- Disable the Secure Boot
- Save your changes and exit

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

- Go to Settings > About
- Your NVIDIA Graphics Card should appear below the Processor field.

## 4) Install CUDA libraries

### Download the CUDA runfile for your operating system

```sh
$ cd ~; mkdir CUDA-Install; cd CUDA-Install
$ wget -nc https://developer.download.nvidia.com/compute/cuda/11.3.0/local_installers/cuda_11.3.0_465.19.01_linux.run
```

### Install dependencies

```sh
$ sudo dnf install gcc-c++ mesa-libGLU-devel libX11-devel libXi-devel libXmu-devel libXcursor-devel libXrandr-devel
$ sudo dnf install freeglut freeglut-devel
```

### Run the CUDA runfile

```sh
$ chmod +x cuda_11.3.0_465.19.01_linux.run
$ sudo ./cuda_11.3.0_465.19.01_linux.run
```

### Configure the installation

- Accept the terms
- Disable the NVIDIA Driver installation
- Enable the installation: CUDA Toolkit, CUDA Samples, CUDA Demo, CUDA Documentation
- Install the package

### Post installation steps

#### Log as root user

```sh
$ su -
```

- If is the first time logging as root, you can set a password to the "root" user with:

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

#### Logout and login again and check the variables

```sh
$ echo $PATH
$ echo $LD_LIBRARY_PATH
```

#### Test the CUDA library

```sh
$ cd /home/<username>/NVIDIA_CUDA-11.1_Samples/1_Utilities/deviceQuery
$ make
$ ./deviceQuery
```

## 5) Configure MonoAlg3D_C

### Download remaining dependencies

```sh
$ sudo dnf install libXcursor-devel libXrandr-devel libXinerama-devel
```

#### If you stop here you already have a working MonoAlg3D binary executable in the /bin folder

### [OPTIONAL] Install and configure MPI for batch simulations

#### Download OpenMPI

```sh
$ cd ~; make OpenMPI-Install; cd OpenMPI-Install
$ wget -nc https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.1.tar.gz
```

#### Build OpenMPI

```sh
$ gunzip -c openmpi-4.1.1.tar.gz | tar xf -
$ cd openmpi-4.1.1
$ ./configure --prefix=/usr/local
$ sudo make all install
```
