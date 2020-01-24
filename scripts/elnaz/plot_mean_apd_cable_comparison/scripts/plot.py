# =================================================================
# Author: Lucas Berg
#
# Program that plot the mean APD
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

def plot_apd_over_a_line (x, apd_1, apd_2, apd_3, apd_3_1, apd_3_2):

    x_numpy = np.asarray(x)
    apd1_numpy = np.asarray(apd_1)
    apd2_numpy = np.asarray(apd_2)
    apd3_numpy = np.asarray(apd_3)
    apd3_1_numpy = np.asarray(apd_3_1)
    apd3_2_numpy = np.asarray(apd_3_2)

    rms_12 = calculate_rms(apd1_numpy,apd2_numpy)
    rms_13 = calculate_rms(apd1_numpy,apd3_numpy)
    rms_13_1 = calculate_rms(apd1_numpy,apd3_1_numpy)
    rms_13_2 = calculate_rms(apd1_numpy,apd3_2_numpy)

    plt.grid()
    plt.plot(x,apd_1,label="SC0",linewidth=1.0)
    plt.plot(x,apd_2,label="SC1.1",linewidth=1.0)
    plt.xlabel("x (um)",fontsize=15)
    plt.ylabel("APD (ms)",fontsize=15)
    plt.title("RMS = %g" % (rms_12),fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.savefig("../outputs/apd-comparison-sc0-sc1_1.png")

    plt.clf()
    plt.grid()
    plt.plot(x,apd_1,label="SC0",linewidth=1.0)
    plt.plot(x,apd_3,label="SC1.2",linewidth=1.0)
    plt.xlabel("x (um)",fontsize=15)
    plt.ylabel("APD (ms)",fontsize=15)
    plt.title("RMS = %g" % (rms_13),fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.savefig("../outputs/apd-comparison-sc0-sc1_2.png")

    plt.clf()
    plt.grid()
    plt.plot(x,apd_1,label="SC0",linewidth=1.0)
    plt.plot(x,apd_3_1,label="SC2.1",linewidth=1.0)
    plt.xlabel("x (um)",fontsize=15)
    plt.ylabel("APD (ms)",fontsize=15)
    plt.title("RMS = %g" % (rms_13_1),fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.savefig("../outputs/apd-comparison-sc0-sc2_1.png")

    plt.clf()
    plt.grid()
    plt.plot(x,apd_1,label="SC0",linewidth=1.0)
    plt.plot(x,apd_3_2,label="SC2.2",linewidth=1.0)
    plt.xlabel("x (um)",fontsize=15)
    plt.ylabel("APD (ms)",fontsize=15)
    plt.title("RMS = %g" % (rms_13_2),fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.savefig("../outputs/apd-comparison-sc0-sc2_2.png")

def main():

    if len(sys.argv) != 6:
        print("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
        print("Usage:> python %s <cells_apd_filename_1> <cells_apd_filename_2> <cells_apd_filename_3>" % sys.argv[0])
        print("                  <cells_apd_filename_4> <cells_apd_filename_5>")
        print("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
        print("<cells_apd_filename> = Input file with the APD of each cell")
        print("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
        print("Example:> python %s ../outputs/mean_apd_sc0.txt ../outputs/mean_apd_sc1_1.txt ../outputs/mean_apd_sc1_2.txt ../outputs/mean_apd_sc2_1.txt ../outputs/mean_apd_sc2_2.txt" % sys.argv[0])
        print("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
        return 1

    cells_apd_filename_1 = sys.argv[1]      # SC0
    cells_apd_filename_2 = sys.argv[2]      # SC1_1
    cells_apd_filename_3 = sys.argv[3]      # SC1_2
    cells_apd_filename_3_1 = sys.argv[4]    # SC2_1
    cells_apd_filename_3_2 = sys.argv[5]    # SC2_2

    # Read the input files as Numpy arrays
    print("Reading inputs files ...") 
    cells_apd_1 = np.genfromtxt(cells_apd_filename_1)                           # The cells APD are unordered !!!    
    cells_apd_2 = np.genfromtxt(cells_apd_filename_2)                           # The cells APD are unordered !!!
    cells_apd_3 = np.genfromtxt(cells_apd_filename_3)                           # The cells APD are unordered !!!
    cells_apd_3_1 = np.genfromtxt(cells_apd_filename_3_1)                       # The cells APD are unordered !!!
    cells_apd_3_2 = np.genfromtxt(cells_apd_filename_3_2)                       # The cells APD are unordered !!!

    x = cells_apd_1[:,0]
    apd_1 = cells_apd_1[:,1]
    apd_2 = cells_apd_2[:,1]
    apd_3 = cells_apd_3[:,1]
    apd_3_1 = cells_apd_3_1[:,1]
    apd_3_2 = cells_apd_3_2[:,1]

    plot_apd_over_a_line(x,apd_1,apd_2,apd_3,apd_3_1,apd_3_2)

if __name__ == "__main__":
	main()
