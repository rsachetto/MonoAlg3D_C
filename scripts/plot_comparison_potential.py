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

    timesteps = timesteps[0:50]
    vms = vms[0:50]

    return timesteps, vms


def plot_transmembrane_potential(t1, v1, t2, v2):

    plt.plot(t1, v1, label="Tissue", c="red", linewidth=2.0)
    plt.plot(t2, v2, label="Purkinje", c="blue", linewidth=2.0)

    plt.grid()
    plt.xlabel("t (ms)",fontsize=15)
    plt.ylabel("V (mV)",fontsize=15)
    #plt.title("Action potential - Cells not within a PMJ site",fontsize=14)
    #plt.title("Action potential - PMJ delay",fontsize=14)
    plt.title("Action potential - PMJ block",fontsize=14)
    plt.legend(loc=0,fontsize=14)
    #plt.savefig("pmj-delay.png")
    plt.savefig("pmj-block.png")
    #plt.savefig("normal-cells.png")
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
