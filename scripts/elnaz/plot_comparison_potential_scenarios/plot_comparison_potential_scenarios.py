# ----------------------------------------------------------------------------------------------------
# Author: Lucas Berg
# Script to plot a comparison between each of the scenarios implemented by Elnaz
# 
# How to use this script:
#
# 1) After running all the simulations and having the results stored in the "output" folder, just run the script "get_aps_for_each_scenario.sh" to captured the Action Potential from a cell over the grid to make the comparison
#
# 2) Run this python script just like the example: 
#	$ python plot_comparison_potential_scenarios.py aps/cell-6240-sc0.txt aps/cell-6240-sc1.txt aps/cell-6240-sc2.txt aps/cell-6240-sc2_1.txt aps/cell-6240-sc2_2.txt aps/cell-6240-sc3.txt aps/cell-6240-sc3_1.txt aps/cell-6240-sc3_2.txt aps/cell-6240-sc3_3.txt 0.02 50
# ----------------------------------------------------------------------------------------------------

import sys
import numpy as np
import matplotlib.pyplot as plt

def calculate_rms (y_true, y_pred):
    loss = np.sqrt(np.mean(np.square(((y_true - y_pred) / y_true))))

    return loss

def read_transmembrane_potential(input_file, dt, print_rate):
    data = np.genfromtxt(input_file)
    n = len(data)
    ms_each_step = dt*print_rate

    end_simulation = (n-1) * ms_each_step
    #print("n = %d" % n)
    #print("ms_each_step = %g" % ms_each_step)
    #print("end_simulation = %g" % end_simulation)

    timesteps = np.arange(0,n)*ms_each_step
    vms = data

    #print(timesteps)
    #print(vms)
    #sys.exit(1)

    return timesteps, vms

def plot_transmembrane_potential_v2 (t1, v1, t2, v2, name1, name2, colorname1, colorname2, cellname):
    for i in range(20):
	print("Working in AP number %d ..." % (i))

	start = i*1000
	end = start + 1000
	rms_value = calculate_rms(v1[start:end],v2[start:end])
	#print("AP %d -- RMS = %g" % (i,rms_value))

	plt.clf()
    	plt.plot(t1, v1, label=name1, c=colorname1, linewidth=1.0)
    	plt.plot(t2, v2, label=name2, c=colorname2, linewidth=0.5)

    	#plt.grid()
    	plt.xlabel("t (ms)",fontsize=15)
    	plt.ylabel("V (mV)",fontsize=15)
    	plt.xlim([start,end])
    	plt.title("Action potential - Cell = %s - RMS = %g" % (cellname,rms_value))
    	plt.legend(loc=1,fontsize=12)
    	#plt.savefig("output/ap-%d-comparison_%s-%s.pdf" % (i,name1,name2))
    	plt.savefig("output/ap-%d-comparison_%s-%s.png" % (i,name1,name2), dpi=300)

def main():

    if len(sys.argv) != 8:
        print("---------------------------------------------------------------------------------------------------------------------")
        print("Usage:> python %s <input_file_1> <input_file_2> <input_file_3> <input_file_4> <input_file_5> <dt> <print_rate>" % sys.argv[0])
        print("---------------------------------------------------------------------------------------------------------------------")
        print("<input_file_1> = Input file with the AP from scenario 0")
        print("<input_file_2> = Input file with the AP from scenario 1.1")
	print("<input_file_3> = Input file with the AP from scenario 1.2")
	print("<input_file_4> = Input file with the AP from scenario 2.1")
	print("<input_file_5> = Input file with the AP from scenario 2.2")
        print("<dt> = Timestep value used for the simulation")
        print("<print_rate> = Print rate used for the simulation")
        print("---------------------------------------------------------------------------------------------------------------------")
        print("Example:> python plot_comparison_potential_scenarios.py aps/cell-6240-sc0.txt aps/cell-6240-sc1_1.txt aps/cell-6240-sc1_2.txt aps/cell-6240-sc2_1.txt aps/cell-6240-sc2_2.txt 0.02 50")
	print("---------------------------------------------------------------------------------------------------------------------")
        return 1

    input_file_1 = sys.argv[1]
    input_file_2 = sys.argv[2]
    input_file_3 = sys.argv[3]
    input_file_4 = sys.argv[4]
    input_file_5 = sys.argv[5]
    dt = float(sys.argv[6])
    print_rate = int(sys.argv[7])

    t1, vm1 = read_transmembrane_potential(input_file_1,dt,print_rate)
    t2, vm2 = read_transmembrane_potential(input_file_2,dt,print_rate)
    t3, vm3 = read_transmembrane_potential(input_file_3,dt,print_rate)
    t4, vm4 = read_transmembrane_potential(input_file_4,dt,print_rate)
    t5, vm5 = read_transmembrane_potential(input_file_5,dt,print_rate)

    #plot_transmembrane_potential(t1,vm1,t2,vm2,t3,vm3,t4,vm4)
    plot_transmembrane_potential_v2(t1,vm1,t2,vm2,"sc0","sc1.1","blue","red","80")
    plot_transmembrane_potential_v2(t1,vm1,t3,vm3,"sc0","sc1.2","blue","red","80")
    plot_transmembrane_potential_v2(t1,vm1,t4,vm4,"sc0","sc2.1","blue","red","80")
    plot_transmembrane_potential_v2(t1,vm1,t5,vm5,"sc0","sc2.2","blue","red","80")


if __name__ == "__main__":
    main()
