import sys
import numpy as np
import matplotlib.pyplot as plt


def read_transmembrane_potential(input_file, dt, print_rate, period):
    data = np.genfromtxt(input_file)
    n = len(data)
    ms_each_step = dt*print_rate

    end_simulation = n / ms_each_step
    end_period = period / ms_each_step

    timesteps = np.arange(0,end_period)*ms_each_step
    vms = data

    num_aps = int( len(vms) / len(timesteps) )

    return timesteps, vms, num_aps, int(ms_each_step)


def plot_transmembrane_potential(t, v, period, num_aps, ms_each_step):

    for i in range(num_aps):
        start = i*period / ms_each_step
        end = start + period / ms_each_step

        #print("Start = %d" % start)
        #print("End = %d" % end)
        
        v_ap = v[start:end]

        plt.plot(t, v_ap, label="AP-%d" % (i+1), linewidth=3.0)

    plt.grid()
    plt.xlabel("t (ms)",fontsize=15)
    plt.ylabel("V (mV)",fontsize=15)
    plt.title("Action potential",fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.savefig("output/multiple-aps.pdf")
    #plt.show()


def main():
	
    if len(sys.argv) != 5:
        print("-------------------------------------------------------------------------")
        print("Usage:> python %s <input_file> <dt> <print_rate> <period>" % sys.argv[0])
        print("-------------------------------------------------------------------------")
        print("<input_file> = Input file with the AP's from each timestep")
        print("<dt> = Timestep value used for the simulation")
        print("<print_rate> = Print rate used for the simulation")
        print("<period> = Basic cycle length for stimulation")
        print("-------------------------------------------------------------------------")
        return 1

    input_file = sys.argv[1]
    dt = float(sys.argv[2])
    print_rate = int(sys.argv[3])
    period = int(sys.argv[4])

    t, vm, num_aps, ms_each_step = read_transmembrane_potential(input_file,dt,print_rate,period)

    #print("Length t = %d" % len(t))
    #print("Length vm = %d" % len(vm))
    #print("num_aps = %d" % num_aps)

    plot_transmembrane_potential(t,vm,period,num_aps,ms_each_step)


if __name__ == "__main__":
    main()
