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
	
def plot_ecg_readings (t1, currents1, t2, currents2, nleads):
	lead_names = ["LA", "RA", "LL", "RL", "V1", "V2", "V3", "V4", "V5", "V6"]
	for i in range(nleads):
		#plt.grid()
		plt.plot(t1, currents1[:,i], label="baseline", c="blue", linewidth=1.0)
		plt.plot(t2, currents2[:,i], label="chronic-diabetes", c="red", linewidth=1.0)
		plt.xlabel("t (ms)",fontsize=15)
		plt.ylabel("Current (mA)",fontsize=15)
		plt.title("ECG reading - Lead %s" % (lead_names[i]),fontsize=14)
		plt.legend(loc=0,fontsize=14)
		#plt.show()
		plt.tight_layout()
		plt.savefig("ecg_lead_%s.png" % (lead_names[i]), dpi=300)
		plt.clf()

def main():
	
	if len(sys.argv) != 3:
		print("-------------------------------------------------------------------------")
		print("Usage:> python %s <input_file_1> <input_file_2>" % sys.argv[0])
		print("-------------------------------------------------------------------------")
		print("<input_file_1> = Input file with the ECG reading from each timestep")
		print("<input_file_2> = Input file with the ECG reading from each timestep")
		print("-------------------------------------------------------------------------")
		return 1

	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]

	t1, currents1, num_leads1 = read_ecg_readings(input_file_1)
	t2, currents2, num_leads2 = read_ecg_readings(input_file_2)

	plot_ecg_readings(t1,currents1,t2,currents2,num_leads1)

if __name__ == "__main__":
	main()
