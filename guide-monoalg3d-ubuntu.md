### GUIDE TO CONFIGURE THE MONOALG_3D ON A UBUNTU 18.04 LTS MACHINE

#### Step 1) Update your packages

```sh
$ sudo apt-get udpate
Reboot your system
```
#### Step 2) Install NVIDIA drivers and essential libraries for CUDA	
    Go to the Software Ubuntu Center and open the advance settings by selecting Software and Updates next to Activities on the top corner of your screen.
	Select Additional Drivers
	Mark the proprietary NVIDIA driver that has been tested and click on Aplly Changes
	Reboot your machine

	Verify if the configuration was correct by going to Settings > Details. Your NVIDIA graphics card name should be on the Graphics section now.

```sh
Install the CUDA libraries by running:
$ sudo apt-get install nvidia-cuda-dev nvidia-cuda-toolkit 
Check if everything is working fine:
$ nvidia-smi
```
#### Step 3) Install other dependencies (optional, for GUI building)
```sh
$ sudo apt-get install libx11-dev libxcursor-dev libxrandr-dev libxinerama-dev libz-dev libxi-dev libglu1-mesa-dev libglvnd-dev
```
