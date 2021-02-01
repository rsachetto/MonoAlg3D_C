import sys
import numpy as np
import matplotlib.pyplot as plt


def read_transmembrane_potential(input_file, dt, print_rate):
    data = np.genfromtxt(input_file)
    n = len(data)
    ms_each_step = dt*print_rate

    end_simulation = n / ms_each_step

    timesteps = np.arange(0,n)*ms_each_step
    vms = data

    return timesteps, vms


def plot_transmembrane_potential(t1, v1, t2, v2):
    #plt.plot(t1[0:150], v1[0:150], label="Purkinje", c="green", linewidth=2.0)
    #plt.plot(t2[0:150], v2[0:150], label="Tissue", c="red", linewidth=2.0)
    #plt.plot(t1[2500:2650], v1[2500:2650], label="Purkinje", c="green", linewidth=2.0)
    #plt.plot(t2[2500:2650], v2[2500:2650], label="Tissue", c="red", linewidth=2.0)
    plt.plot(t1[0:2000], v1[0:2000], label="Trovato_2020", c="green", linewidth=2.0)
    plt.plot(t2[0:2000], v2[0:2000], label="ToRORd_2019", c="red", linewidth=2.0)
    #plt.plot(t1[2500:4500], v1[2500:4500], label="Purkinje", c="green", linewidth=2.0)
    #plt.plot(t2[2500:4500], v2[2500:4500], label="Tissue", c="red", linewidth=2.0)
    #plt.plot(t1, v1, label="Purkinje", c="green", linewidth=2.0)
    #plt.plot(t2, v2, label="Tissue", c="red", linewidth=2.0, alpha=0.75)

    #plt.grid()
    #plt.xlim([500,600])
    plt.xlabel("t (ms)",fontsize=15)
    plt.ylabel("V (mV)",fontsize=15)
    plt.title("Action potential traces",fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.savefig("trovato_torord_comparison.png")
    #plt.show()


def main():
	
    if len(sys.argv) != 5:
        print("-------------------------------------------------------------------------")
        print("Usage:> python %s <input_file_1> <input_file_2> <dt> <print_rate>" % sys.argv[0])
        print("-------------------------------------------------------------------------")
        print("<input_file_1> = Input file with the AP from the first simulation")
        print("<input_file_2> = Input file with the AP from the second simulation")
        print("<dt> = Timestep value used for the simulation")
        print("<print_rate> = Print rate used for the simulation")
        print("-------------------------------------------------------------------------")
        return 1

    input_file_1 = sys.argv[1]
    input_file_2 = sys.argv[2]
    dt = float(sys.argv[3])
    print_rate = int(sys.argv[4])

    t1, vm1 = read_transmembrane_potential(input_file_1,dt,print_rate)
    t2, vm2 = read_transmembrane_potential(input_file_2,dt,print_rate)

    plot_transmembrane_potential(t1,vm1,t2,vm2)


if __name__ == "__main__":
    main()
