import sys
import numpy as np
import matplotlib.pyplot as plt

def read_state_vector(input_file):
    sv = np.genfromtxt(input_file)

    timesteps = sv[:,0]
    sv = sv[:,1:]

    return timesteps, sv

def plot_all_state_vector(t, sv):
    (nlin, ncol) = sv.shape

    for i in range(ncol):
        plt.plot(t,sv[:,i],label="sv[%d]" % (i))
        plt.legend(loc=0,fontsize=14)
        #plt.show()
        plt.savefig("traces/torord/sv_%d.svg" % (i))
        plt.clf()

    for i in range(ncol):
        print("sv[%d] = %g;" % (i,sv[nlin-1][i]))

def plot_state_vector(t, sv, col):
    plt.plot(t,sv[:,col],label="sv[%d]" % (col))
    plt.legend(loc=0,fontsize=14)
    plt.show()

def main():
	
    if len(sys.argv) != 2:
        print("-------------------------------------------------------------------------")
        print("Usage:> python %s <input_file>" % sys.argv[0])
        print("-------------------------------------------------------------------------")
        print("<input_file> = Input file with the state-vector data")
        print("-------------------------------------------------------------------------")
        return 1

    input_file = sys.argv[1]

    t, sv = read_state_vector(input_file)

    plot_all_state_vector(t,sv)

    # This will plot only the transmembrane potential
    #plot_state_vector(t,sv,0)

if __name__ == "__main__":
	main()
