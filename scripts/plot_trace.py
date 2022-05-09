# ==============================================================================================================================
# Author: Lucas Berg
# ==============================================================================================================================
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def main ():
    if len(sys.argv) != 5:
        print("----------------------------------------------------------------------------------------------------------------------")
        print("Usage:> python %s <input_file> <tmax> <dt> <printrate>" % sys.argv[0])
        print("----------------------------------------------------------------------------------------------------------------------")
        print("<input_file> = Input file with the state-vector data")
        print("<tmax> = Maximum simulation time")
        print("<dt> = Timestep used in the simulation")
        print("<printrate> = Print rate used in the simulation")
        print("----------------------------------------------------------------------------------------------------------------------")
        return 1

    input_file = sys.argv[1]
    tmax = float(sys.argv[2])
    dt = float(sys.argv[3])
    printrate = int(sys.argv[4])

    sv = np.genfromtxt(input_file)

    plt.plot(sv[:,0],sv[:,1])
    plt.show()

if __name__ == "__main__":
	main()