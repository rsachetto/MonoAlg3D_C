import sys
import numpy as np
import matplotlib.pyplot as plt

def plot_transmembrane_potential(t, v):
	plt.plot(t, v, label="Vm", c="black", linewidth=1.0)
	plt.xlabel("t (ms)",fontsize=15)
	plt.ylabel("V (mV)",fontsize=15)
	plt.title("Action potential",fontsize=14)
	plt.legend(loc=0,fontsize=14)
	plt.show()
	#plt.savefig("ap.pdf")

def main():
	
	if len(sys.argv) != 2:
		print("-------------------------------------------------------------------------")
		print("Usage:> python %s <input_file>" % sys.argv[0])
		print("-------------------------------------------------------------------------")
		print("<input_file> = Input file with the state vector from each timestep")
		print("-------------------------------------------------------------------------")
		return 1

	input_file = sys.argv[1]

	sv = np.genfromtxt(input_file)

	t = sv[:,0]
	vm = sv[:,1]

	plot_transmembrane_potential(t,vm)

if __name__ == "__main__":
	main()
