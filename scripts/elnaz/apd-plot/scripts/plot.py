# =================================================================
# Author: Lucas Berg
#
# Program that plot the APD over a line
# This script only works if the grid is a plain tissue
# =================================================================

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import interpolate

def calculate_rms (y_true, y_pred):
    loss = np.sqrt(np.mean(np.square(((y_true - y_pred) / y_true))))

    return loss

def plot_apd_over_a_line (x, apd_1, apd_2, apd_3, apd_4, apd_4_1):

    x_numpy = np.asarray(x)
    apd1_numpy = np.asarray(apd_1)
    apd2_numpy = np.asarray(apd_2)
    apd3_numpy = np.asarray(apd_3)
    apd4_numpy = np.asarray(apd_4)
    apd4_1_numpy = np.asarray(apd_4_1)

    rms_12 = calculate_rms(apd1_numpy,apd2_numpy)
    rms_13 = calculate_rms(apd1_numpy,apd3_numpy)
    rms_14 = calculate_rms(apd1_numpy,apd4_numpy)
    rms_14_1 = calculate_rms(apd1_numpy,apd4_1_numpy)

    plt.grid()
    plt.plot(x,apd_1,label="SC1",linewidth=1.0)
    plt.plot(x,apd_2,label="SC2",linewidth=1.0)
    plt.xlabel("x (um)",fontsize=15)
    plt.ylabel("APD (ms)",fontsize=15)
    plt.title("RMS = %g" % (rms_12),fontsize=14)
    plt.legend(loc=0,fontsize=14)
    #plt.savefig("output/apd-comparison-sc1-sc2.pdf")
    plt.savefig("../outputs/apd-comparison-sc1-sc2.png")

    plt.clf()
    plt.grid()
    plt.plot(x,apd_1,label="SC1",linewidth=1.0)
    plt.plot(x,apd_3,label="SC3",linewidth=1.0)
    plt.xlabel("x (um)",fontsize=15)
    plt.ylabel("APD (ms)",fontsize=15)
    plt.title("RMS = %g" % (rms_13),fontsize=14)
    plt.legend(loc=0,fontsize=14)
    #plt.savefig("output/apd-comparison-sc1-sc3.pdf")
    plt.savefig("../outputs/apd-comparison-sc1-sc3.png")

    plt.clf()
    plt.grid()
    plt.plot(x,apd_1,label="SC1",linewidth=1.0)
    plt.plot(x,apd_4,label="SC4",linewidth=1.0)
    plt.xlabel("x (um)",fontsize=15)
    plt.ylabel("APD (ms)",fontsize=15)
    plt.title("RMS = %g" % (rms_14),fontsize=14)
    plt.legend(loc=0,fontsize=14)
    #plt.savefig("output/apd-comparison-sc1-sc4.pdf")
    plt.savefig("../outputs/apd-comparison-sc1-sc4.png")

    plt.clf()
    plt.grid()
    plt.plot(x,apd_1,label="SC1",linewidth=1.0)
    plt.plot(x,apd_4_1,label="SC4_1",linewidth=1.0)
    plt.xlabel("x (um)",fontsize=15)
    plt.ylabel("APD (ms)",fontsize=15)
    plt.title("RMS = %g" % (rms_14_1),fontsize=14)
    plt.legend(loc=0,fontsize=14)
    #plt.savefig("output/apd-comparison-sc1-sc4_1.pdf")
    plt.savefig("../outputs/apd-comparison-sc1-sc4_1.png")


def main():
	
    if len(sys.argv) != 6:
        print("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
        print("Usage:> python %s <cells_apd_filename_1> <cells_apd_filename_2> <cells_apd_filename_3>" % sys.argv[0])
        print("                  <cells_apd_filename_4> <cells_apd_filename_5>")
        print("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
        print("<cells_apd_filename> = Input file with the APD of each cell")
        print("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
        print("Example:> python %s inputs/" % sys.argv[0])
        return 1

    cells_apd_filename_1 = sys.argv[1]
    cells_apd_filename_2 = sys.argv[2]
    cells_apd_filename_3 = sys.argv[3]
    cells_apd_filename_4 = sys.argv[4]
    cells_apd_filename_4_1 = sys.argv[5]

    # Read the input files as Numpy arrays
    print("Reading inputs files ...") 
    cells_apd_1 = np.genfromtxt(cells_apd_filename_1)                       # The cells APD are unordered !!!
    cells_apd_2 = np.genfromtxt(cells_apd_filename_2)                       # The cells APD are unordered !!!
    cells_apd_3 = np.genfromtxt(cells_apd_filename_3)                       # The cells APD are unordered !!!
    cells_apd_4 = np.genfromtxt(cells_apd_filename_4)                       # The cells APD are unordered !!!
    cells_apd_5 = np.genfromtxt(cells_apd_filename_4_1)                       # The cells APD are unordered !!!

    x = cells_apd_1[:,0]
    apd_1 = cells_apd_1[:,1]
    apd_2 = cells_apd_2[:,1]
    apd_3 = cells_apd_3[:,1]
    apd_4 = cells_apd_4[:,1]
    apd_4_1 = cells_apd_5[:,1]
    

    plot_apd_over_a_line(x,apd_1,apd_2,apd_3,apd_4,apd_4_1)
    

if __name__ == "__main__":
	main()
