## GUIDE TO CONFIGURE MONOALG_3D ON A WINDOWS 11 MACHINE

----------------------------------------------------------

This guide provides instructions to configure the MonoAlg3D solver in a fresh instalation of Windows 11.

## 1) Enable WSL2

- Enable Virtualisation using the Boot Menu
- Open the "Enable and Disable Windows resources"
- Activate "Linux Subsystem for Windows" and "Virtual Machine Platform"
- Reboot your machine
- Open a PowerShell with Administrator privileges
- Execute the command:
```sh
	$ wsl.exe --install
```
- Create a username for the Ubuntu virtual machine
- Perform a system update inside the Ubuntu machine
```sh
	$ sudo apt-get update
```

## 2) Install dependencies

In your Ubuntu machine do the following:

- Install GCC/G++ compilers
```sh
	$ sudo apt-get install build-essential
```
- Install other packages
```sh
	$ sudo apt-get install libx11-dev libxcursor-dev libxrandr-dev libxinerama-dev libz-dev libxi-dev libglu1-mesa-dev libglvnd-dev
```
- Check if you have the NVIDIA driver configured and enabled
```sh
	$ nvidia-smi
```

## 3) Install CUDA

```sh
	$ cd ~; mkdir CUDA-Install; cd CUDA-Install
	$ wget -nc https://developer.download.nvidia.com/compute/cuda/12.6.3/local_installers/cuda_12.6.3_560.35.05_linux.run
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

#### Set the PATH and LD_LIBRARY_PATH for CUDA

```sh
$ cd ~
$ nano .bashrc
```

- Open your _"~/.bashrc"_ file and include the following lines:

```sh
PATH=$PATH:/usr/local/cuda/include:/usr/local/cuda/bin
LD_LIBRARY_PATH=/usr/local/lib:/usr/local/cuda/lib64
export PATH
export LD_LIBRARY_PATH
```

#### Logout and login again or execute _"source .bashrc"_. Check the variables

```sh
$ echo $PATH
$ echo $LD_LIBRARY_PATH
```

- The new paths must appear.

#### Test the CUDA library

```sh
$ git clone https://github.com/bergolho1337/CUDA.git
$ cd CUDA/DeviceInfo
$ make
$ ./deviceInfo
```

- This should print some information about your GPU hardware.

## 4) Configure MonoAlg3D_C

### Download remaining dependencies

```sh
$ sudo apt-get install libx11-dev libxcursor-dev libxrandr-dev libxinerama-dev libz-dev libxi-dev libglu1-mesa-dev libglvnd-dev
```

### Update the _"bsbash/find_functions.sh"_ file with the following:

```sh
CUDA_LIBRARY_PATH="/usr/local/cuda/lib64"
CUDA_MATH_LIBRARY_PATH=""
CUDA_INCLUDE_PATH="/usr/local/cuda/include"
NVCC=""
CUDA_FOUND=""
```

### Build and compile MonoAlg3D

```sh
$ ./build.sh -f simulator
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
