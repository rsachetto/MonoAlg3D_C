import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_data (data):
    plt.grid()
    plt.plot(data[:,0], data[:,1], label="PMJ Delay tuning - NminPMJ=60", c="black", linewidth=2.0, marker='o')
    plt.xlabel(r"$R_{PMJ}$ ($M \Omega$)",fontsize=15)
    plt.ylabel("PMJ delay (ms)",fontsize=15)
    plt.ylim([0,20])
    plt.title("Trovato_TT3 - PMJ delay",fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.show()
    #plt.savefig("pmj_delay_tt3_trovato.pdf")

def main():
	
    if len(sys.argv) != 2:
        print("-------------------------------------------------------------------------")
        print("Usage:> python %s <input_file>" % sys.argv[0])
        print("-------------------------------------------------------------------------")
        print("<input_file> = Input file with measurements")
        print("-------------------------------------------------------------------------")
        print("Example:> python plot.py ../rpmj_delay.txt")
        print("-------------------------------------------------------------------------")
        return 1
    
    input_file = sys.argv[1]
    data = np.genfromtxt(input_file)
    plot_data(data)

if __name__ == "__main__":
	main()
