'''
====================================================================================================================================
Version: 12/01/2024
Author: Lucas Berg (@bergolho)
====================================================================================================================================
Script to read and plot the ECG data coming from a MonoAlg3D simulation using Matplotlib.
------------------------------------------------------------------------------------------------------------------------------------
How to use:> python plot_ecg.py <path_to_output_simulation_folder>/ecg.txt
====================================================================================================================================
'''

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

	fig, axs = plt.subplots(nleads, 1, sharex=True, sharey=False)
	fig.text(0.5, 0.04, 'Time (ms)', ha='center', fontsize=10)
	fig.text(0.02, 0.5, 'Current (mA)', va='center', rotation='vertical', fontsize=10)
	for i in range(nleads):
		axs[i].plot(t, currents[:,i], label="lead_%d" % (i), c="black", linewidth=2.0)
		axs[i].set_title("ECG reading - Lead %d" % (i),fontsize=14)
	
	#plt.show()
	#plt.savefig("ecg.pdf")
	plt.savefig("ecg.png", dpi=300)
	print("[+] ECG readings saved in the current folder with the name 'ecg.png'!")

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
