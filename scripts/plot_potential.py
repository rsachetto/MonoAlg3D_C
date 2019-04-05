import sys
import numpy as np
import matplotlib.pyplot as plt


def read_transmembrane_potential(input_file, dt, print_rate):
	data = np.genfromtxt(input_file)
	n = len(data)

	timesteps = np.arange(0,n)*dt*print_rate
	vms = data
	
	return timesteps, vms


def plot_transmembrane_potential(t, v):
	plt.grid()
	plt.plot(t, v, label="Vm", c="black", linewidth=3.0)
	plt.xlabel("t (ms)",fontsize=15)
	plt.ylabel("V (mV)",fontsize=15)
	plt.title("Action potential",fontsize=14)
	plt.legend(loc=2,fontsize=14)
	#plt.show()
	plt.savefig("ap.pdf")


def main():
	
	if len(sys.argv) != 4:
		print("-------------------------------------------------------------------------")
		print("Usage:> python %s <input_file> <dt> <print_rate>" % sys.argv[0])
		print("-------------------------------------------------------------------------")
		print("<input_file> = Input file with the AP's from each timestep")
		print("<dt> = Timestep value used for the simulation")
		print("<print_rate> = Print rate used for the simulation")
		print("-------------------------------------------------------------------------")
		return 1

	input_file = sys.argv[1]
	dt = float(sys.argv[2])
	print_rate = int(sys.argv[3])

	t, vm = read_transmembrane_potential(input_file,dt,print_rate)

	plot_transmembrane_potential(t,vm)


if __name__ == "__main__":
	main()
