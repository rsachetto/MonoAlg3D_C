# ==============================================================================================================================
# Author: Lucas Berg
# ==============================================================================================================================
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def main ():
    if len(sys.argv) != 6:
        print("----------------------------------------------------------------------------------------------------------------------")
        print("Usage:> python %s <input_file_1> <input_file_2> <tmax> <dt> <printrate>" % sys.argv[0])
        print("----------------------------------------------------------------------------------------------------------------------")
        print("<input_file_1> = Input file with the state-vector data 1")
        print("<input_file_2> = Input file with the state-vector data 2")
        print("<tmax> = Maximum simulation time")
        print("<dt> = Timestep used in the simulation")
        print("<printrate> = Print rate used in the simulation")
        print("----------------------------------------------------------------------------------------------------------------------")
        return 1

    input_file_1 = sys.argv[1]
    input_file_2 = sys.argv[2]
    tmax = float(sys.argv[3])
    dt = float(sys.argv[4])
    printrate = int(sys.argv[5])

    sv_1 = np.genfromtxt(input_file_1)
    sv_2 = np.genfromtxt(input_file_2)

    vm_1 = sv_1[:,1]
    vm_2 = sv_2[:,1]

    #plt.grid()
    plt.plot(sv_1[:,0], sv_1[:,1], label="pk", c="red", linewidth=2.0)
    plt.plot(sv_2[:,0], sv_2[:,1], label="myo", c="blue", linewidth=2.0)
    plt.xlabel("Time (ms)",fontsize=15)
    plt.ylabel("Transmembrane potential (mV)",fontsize=15)
    plt.title("Action potential",fontsize=14)
    plt.legend(loc=0,fontsize=14)
    #plt.show()
    plt.savefig("purkinje_coupling_aps.pdf")

if __name__ == "__main__":
	main()
