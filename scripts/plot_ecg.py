import sys
import numpy as np
import matplotlib.pyplot as plt

def read_ecg_readings (input_file):
	data = np.genfromtxt(input_file)
	nlin, ncol = np.shape(data)
	timesteps = data[:,0]
	currents = data[:,1:]
	num_leads = ncol-1
	
	return timesteps, currents, num_leads
	
def plot_ecg_readings (t, currents, nleads):
	for i in range(nleads):
		#plt.grid()
		plt.plot(t, currents[:,i], label="lead_%d" % (i), c="black", linewidth=3.0)
		plt.xlabel("t (ms)",fontsize=15)
		plt.ylabel("Current (mA)",fontsize=15)
		plt.title("ECG reading - Lead %d" % (i),fontsize=14)
		plt.legend(loc=0,fontsize=14)
		plt.show()
		#plt.savefig("ecg_lead_%d.pdf" % (i))

def main():
	
	if len(sys.argv) != 2:
		print("-------------------------------------------------------------------------")
		print("Usage:> python %s <input_file>" % sys.argv[0])
		print("-------------------------------------------------------------------------")
		print("<input_file> = Input file with the ECG reading from each timestep")
		print("-------------------------------------------------------------------------")
		return 1

	input_file = sys.argv[1]

	t, currents, num_leads = read_ecg_readings(input_file)

	plot_ecg_readings(t,currents,num_leads)

if __name__ == "__main__":
	main()
