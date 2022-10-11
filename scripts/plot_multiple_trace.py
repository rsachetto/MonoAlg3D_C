# ==============================================================================================================================
# Author: Lucas Berg
# ==============================================================================================================================
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def main ():
    if len(sys.argv) != 7:
        print("----------------------------------------------------------------------------------------------------------------------")
        print("Usage:> python %s <input_file_1> <input_file_2> <input_file_3> <tmax> <dt> <printrate>" % sys.argv[0])
        print("----------------------------------------------------------------------------------------------------------------------")
        print("<input_file_1> = Input file with the state-vector data 1")
        print("<input_file_2> = Input file with the state-vector data 2")
        print("<input_file_3> = Input file with the state-vector data 3")
        print("<tmax> = Maximum simulation time")
        print("<dt> = Timestep used in the simulation")
        print("<printrate> = Print rate used in the simulation")
        print("----------------------------------------------------------------------------------------------------------------------")
        return 1

    input_file_1 = sys.argv[1]
    input_file_2 = sys.argv[2]
    input_file_3 = sys.argv[3]
    tmax = float(sys.argv[4])
    dt = float(sys.argv[5])
    printrate = int(sys.argv[6])

    sv_1 = np.genfromtxt(input_file_1)
    sv_2 = np.genfromtxt(input_file_2)
    sv_3 = np.genfromtxt(input_file_3)

    vm_1 = sv_1[0:25000,1]
    vm_2 = sv_2[0:25000,1]
    vm_3 = sv_3[0:25000,1]

    t = []
    num_timesteps = len(sv_1[:,0])
    for i in range(25000):
        t.append(dt*i)
    t = np.array(t)
    
    #plt.grid()
    plt.plot(t,vm_1, label="endocardial", c="blue", linewidth=2.0)
    plt.plot(t,vm_2, label="epicardial", c="red", linewidth=2.0)
    plt.plot(t,vm_3, label="midmyocardial", c="goldenrod", linewidth=2.0)
    plt.xlabel("Time (ms)",fontsize=15)
    plt.ylabel("Transmembrane potential (mV)",fontsize=15)
    plt.title("Action potential - ToRORd_dynCl - Transmural (2020)",fontsize=14)
    plt.legend(loc=0,fontsize=14)
    #plt.show()
    plt.savefig("action_potential_ToRORd_dynCl_2020_transmural_bcl_1000ms.pdf")

if __name__ == "__main__":
	main()