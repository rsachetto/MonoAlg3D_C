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

def plot_ecg_readings_dti003(t, currents, nleads):

	x_arr, y_arr = [], []
	for k in range(nleads):
		x, y = [], []
		for i in range(len(t)):
			x.append(t[i])
			y.append(float(currents[i][k]))
		x_arr.append(x)
		y_arr.append(y)

	fig, ((ax1, ax2, ax3, ax4, ax5), (ax6, ax7, ax8, ax9, ax10)) = plt.subplots(2, 5, sharex=True, sharey=True)
	fig.suptitle('Calculated ECG - DTI003', fontsize=14)
	fig.text(0.5, 0.04, 'Time (ms)', ha='center', fontsize=10)
	fig.text(0.04, 0.5, 'Measured value', va='center', rotation='vertical', fontsize=10)
	ax1.plot(x_arr[0],y_arr[0],c="black",linewidth=1.0)
	ax1.set_title("Electrode 1")
	ax2.plot(x_arr[1],y_arr[1],c="black",linewidth=1.0)
	ax2.set_title("Electrode 2")
	ax3.plot(x_arr[2],y_arr[2],c="black",linewidth=1.0)
	ax3.set_title("Electrode 3")
	ax4.plot(x_arr[3],y_arr[3],c="black",linewidth=1.0)
	ax4.set_title("Electrode 4")
	ax5.plot(x_arr[4],y_arr[4],c="black",linewidth=1.0)
	ax5.set_title("Electrode 5")
	ax6.plot(x_arr[5],y_arr[5],c="black",linewidth=1.0)
	ax6.set_title("Electrode 6")
	ax7.plot(x_arr[6],y_arr[6],c="black",linewidth=1.0)    
	ax7.set_title("Electrode 7")
	ax8.plot(x_arr[7],y_arr[7],c="black",linewidth=1.0)
	ax8.set_title("Electrode 8")
	ax9.plot(x_arr[8],y_arr[8],c="black",linewidth=1.0)
	ax9.set_title("Electrode 9")
	ax10.plot(x_arr[9],y_arr[9],c="black",linewidth=1.0)
	ax10.set_title("Electrode 10")

	for ax in fig.get_axes():
		ax.label_outer()

	#plt.show()
	plt.savefig("ecg_readings_dti003.png",dpi=150)
	
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

	#plot_ecg_readings(t,currents,num_leads)
	plot_ecg_readings_dti003(t,currents,num_leads)

if __name__ == "__main__":
	main()
