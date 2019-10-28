# Author: Lucas Berg and Elnaz
# Script to plot a comparison between each of the 4 scenarios from Elnaz`s experiments


import sys
import numpy as np
import matplotlib.pyplot as plt


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


def plot_transmembrane_potential(t1, v1, t2, v2, t3, v3, t4, v4):

    plt.plot(t1, v1, label="SC0", c="blue", linewidth=1.0)
    plt.plot(t2, v2, label="SC1", c="red", linewidth=1.0)
    plt.plot(t3, v3, label="SC2", c="green", linewidth=1.0)
    plt.plot(t4, v4, label="SC3", c="purple", linewidth=1.0)

    plt.grid()
    plt.xlim(1200.0,1800.0)
    plt.xlabel("t (ms)",fontsize=15)
    plt.ylabel("V (mV)",fontsize=15)
    plt.title("Action potential - Cell = 5876",fontsize=14)
    plt.legend(loc=1,fontsize=12)
    plt.savefig("output/comparison-aps.pdf")
    #plt.show()


def main():
	
    if len(sys.argv) != 7:
        print("---------------------------------------------------------------------------------------------------------------------")
        print("Usage:> python %s <input_file_1> <input_file_2> <input_file_3> <input_file_4> <dt> <print_rate>" % sys.argv[0])
        print("---------------------------------------------------------------------------------------------------------------------")
        print("<input_file_1> = Input file with the AP from scenario 0")
        print("<input_file_2> = Input file with the AP from scenario 1")
	print("<input_file_3> = Input file with the AP from scenario 2")
 	print("<input_file_4> = Input file with the AP from scenario 3")
        print("<dt> = Timestep value used for the simulation")
        print("<print_rate> = Print rate used for the simulation")
        print("---------------------------------------------------------------------------------------------------------------------")
        return 1

    input_file_1 = sys.argv[1]
    input_file_2 = sys.argv[2]
    input_file_3 = sys.argv[3]
    input_file_4 = sys.argv[4]
    dt = float(sys.argv[5])
    print_rate = int(sys.argv[6])

    t1, vm1 = read_transmembrane_potential(input_file_1,dt,print_rate)
    t2, vm2 = read_transmembrane_potential(input_file_2,dt,print_rate)
    t3, vm3 = read_transmembrane_potential(input_file_3,dt,print_rate)
    t4, vm4 = read_transmembrane_potential(input_file_4,dt,print_rate)

    plot_transmembrane_potential(t1,vm1,t2,vm2,t3,vm3,t4,vm4)


if __name__ == "__main__":
    main()
