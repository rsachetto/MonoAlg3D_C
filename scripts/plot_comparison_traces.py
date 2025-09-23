# ==============================================================================================================================
# Author: Lucas Berg
# ==============================================================================================================================
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def main ():
    if len(sys.argv) != 3:
        print("----------------------------------------------------------------------------------------------------------------------")
        print("Usage:> python %s <input_file_1> <input_file_2>" % sys.argv[0])
        print("----------------------------------------------------------------------------------------------------------------------")
        print("<input_file_1> = Input file with the state-vector data from the reference")
        print("<input_file_2> = Input file with the state-vector data from the approxition")
        print("----------------------------------------------------------------------------------------------------------------------")
        return 1

    input_file_1 = sys.argv[1]
    input_file_2 = sys.argv[2]

    sv_1 = np.genfromtxt(input_file_1)
    sv_2 = np.genfromtxt(input_file_2)

    t_1 = sv_1[:,0]
    vm_1 = sv_1[:,1]
    cai_1 = sv_1[:,2]
    #ta_1 = sv_1[:,3]

    t_2 = sv_2[:,0]
    vm_2 = sv_2[:,1]
    cai_2 = sv_2[:,2]
    #ta_2 = sv_2[:,3]

    plt.plot(t_1,vm_1, label="Matlab", c="red", linewidth=1.0)
    plt.plot(t_2-249000,vm_2, label="MonoAlg3D", c="black", linewidth=0.5)
    #plt.plot(t_1[1130:],vm_1[1130:], label="Matlab", c="red", linewidth=1.0)
    #plt.plot(t_2[600000:],vm_2[600000:], label="MonoAlg3D", c="black", linewidth=0.5)    # Last of 4 beats
    #plt.plot(t_1,vm_1, label="Matlab", c="red", linewidth=1.0)  # Last of 20 beats
    #plt.plot(t_2[3800000:]-19000,vm_2[3800000:], label="MonoAlg3D", c="black", linewidth=0.5)  # Last of 20 beats
    plt.xlabel("Time (ms)",fontsize=15)
    plt.ylabel("Transmembrane potential (mV)",fontsize=15)
    #plt.title("Vm - Euler adaptive [reltol=1e-9]",fontsize=14)
    plt.title("Vm - Rush-Larsen [dt=0.005ms] - Beat 250",fontsize=14)
    plt.legend(loc=0,fontsize=14)
    #plt.show()
    plt.tight_layout()
    plt.savefig("action_potential.png", dpi=300)

    plt.clf()
    plt.plot(t_1,cai_1, label="Matlab", c="red", linewidth=1.0)
    plt.plot(t_2-249000,cai_2, label="MonoAlg3D", c="black", linewidth=0.5)
    #plt.plot(t_1[1130:],cai_1[1130:], label="Matlab", c="red", linewidth=1.0)
    #plt.plot(t_2[600000:],cai_2[600000:], label="MonoAlg3D", c="black", linewidth=0.5)  # Last of 4 beats
    #plt.plot(t_1,cai_1, label="Matlab", c="red", linewidth=1.0)   # Last of 20 beats
    #plt.plot(t_2[3800000:]-19000,cai_2[3800000:], label="MonoAlg3D", c="black", linewidth=0.5)  # Last of 20 beats
    plt.xlabel("Time (ms)",fontsize=15)
    plt.ylabel("[mM]",fontsize=15)
    #plt.title("Calcium - Euler adaptive [reltol=1e-9]",fontsize=14)
    plt.title("Calcium - Rush-Larsen [dt=0.005ms] - Beat 250",fontsize=14)
    plt.legend(loc=0,fontsize=14)
    #plt.show()
    plt.tight_layout()
    plt.savefig("cai.png", dpi=300)

    #plt.clf()
    #plt.plot(t_1,ta_1, label="ref", c="red", linewidth=1.0)
    #plt.plot(t_2,ta_2, label="aprox", c="black", linewidth=0.5)
    #plt.xlabel("Time (ms)",fontsize=15)
    #plt.ylabel("[mM]",fontsize=15)
    #plt.title("Active tension",fontsize=14)
    #plt.legend(loc=0,fontsize=14)
    #plt.show()
    #plt.savefig("ta.png", dpi=300)


if __name__ == "__main__":
	main()
