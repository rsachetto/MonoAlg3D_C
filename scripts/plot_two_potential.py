import sys
import numpy as np
import matplotlib.pyplot as plt


def read_transmembrane_potential(input_file, dt, print_rate):
	data = np.genfromtxt(input_file)
	n = len(data)

	timesteps = np.arange(0,n)*dt*print_rate
	vms = data

	timesteps = timesteps[500:]
	vms = vms[500:]

	return timesteps, vms


def plot_transmembrane_potential(t, v1, v2):
	#plt.grid()
	plt.plot(t, v1, label="Tissue", c="blue", linewidth=3.0)
	plt.plot(t, v2, label="Purkinje", c="orange", linewidth=3.0)
	plt.xlabel("t (ms)",fontsize=15)
	plt.ylabel("V (mV)",fontsize=15)
	plt.title("Action potential",fontsize=14)
	plt.legend(loc=0,fontsize=14)
	#plt.show()
	#plt.savefig("two_ap.pdf")
	plt.savefig("two_ap.png",dpi=300)


def main():
	
	if len(sys.argv) != 5:
		print("-------------------------------------------------------------------------")
		print("Usage:> python %s <input_file_1> <input_file_2> <dt> <print_rate>" % sys.argv[0])
		print("-------------------------------------------------------------------------")
		print("<input_file_1> = Input file with the AP's from cell 1")
		print("<input_file_2> = Input file with the AP's from cell 2")
		print("<dt> = Timestep value used for the simulation")
		print("<print_rate> = Print rate used for the simulation")
		print("-------------------------------------------------------------------------")
		return 1

	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]
	dt = float(sys.argv[3])
	print_rate = int(sys.argv[4])

	t1, vm1 = read_transmembrane_potential(input_file_1,dt,print_rate)
	t2, vm2 = read_transmembrane_potential(input_file_2,dt,print_rate)

	plot_transmembrane_potential(t1,vm1,vm2)


if __name__ == "__main__":
	main()
