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
    plt.savefig("output/apd-comparison-sc1-sc2.png")

    plt.clf()
    plt.grid()
    plt.plot(x,apd_1,label="SC1",linewidth=1.0)
    plt.plot(x,apd_3,label="SC3",linewidth=1.0)
    plt.xlabel("x (um)",fontsize=15)
    plt.ylabel("APD (ms)",fontsize=15)
    plt.title("RMS = %g" % (rms_13),fontsize=14)
    plt.legend(loc=0,fontsize=14)
    #plt.savefig("output/apd-comparison-sc1-sc3.pdf")
    plt.savefig("output/apd-comparison-sc1-sc3.png")

    plt.clf()
    plt.grid()
    plt.plot(x,apd_1,label="SC1",linewidth=1.0)
    plt.plot(x,apd_4,label="SC4",linewidth=1.0)
    plt.xlabel("x (um)",fontsize=15)
    plt.ylabel("APD (ms)",fontsize=15)
    plt.title("RMS = %g" % (rms_14),fontsize=14)
    plt.legend(loc=0,fontsize=14)
    #plt.savefig("output/apd-comparison-sc1-sc4.pdf")
    plt.savefig("output/apd-comparison-sc1-sc4.png")

    plt.clf()
    plt.grid()
    plt.plot(x,apd_1,label="SC1",linewidth=1.0)
    plt.plot(x,apd_4_1,label="SC4_1",linewidth=1.0)
    plt.xlabel("x (um)",fontsize=15)
    plt.ylabel("APD (ms)",fontsize=15)
    plt.title("RMS = %g" % (rms_14_1),fontsize=14)
    plt.legend(loc=0,fontsize=14)
    #plt.savefig("output/apd-comparison-sc1-sc4_1.pdf")
    plt.savefig("output/apd-comparison-sc1-sc4_1.png")


def main():
	
    if len(sys.argv) != 7:
        print("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
        print("Usage:> python %s <sorted_cells_position_filename> <cells_apd_filename_1> <cells_apd_filename_2> <cells_apd_filename_3> <cells_apd_filename_4> <cells_apd_filename_5>" % sys.argv[0])
        print("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
        print("<sorted_cells_position_filename> = Input file with the sorted positions of each cell")
        print("<cells_apd_filename> = Input file with the APD of each cell")
        print("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
        print("Example:> python %s inputs/cells_positions_inside_region.txt outputs/cells-apd-inside-region.txt" % sys.argv[0])
        return 1

    sorted_cells_position_filename = sys.argv[1]
    cells_apd_filename_1 = sys.argv[2]
    cells_apd_filename_2 = sys.argv[3]
    cells_apd_filename_3 = sys.argv[4]
    cells_apd_filename_4 = sys.argv[5]
    cells_apd_filename_4_1 = sys.argv[6]

    # Read the input files as Numpy arrays
    print("Reading inputs files ...") 
    sorted_cells_positions = np.genfromtxt(sorted_cells_position_filename)
    cells_apd_1 = np.genfromtxt(cells_apd_filename_1)                       # The cells APD are unordered !!!
    cells_apd_2 = np.genfromtxt(cells_apd_filename_2)                       # The cells APD are unordered !!!
    cells_apd_3 = np.genfromtxt(cells_apd_filename_3)                       # The cells APD are unordered !!!
    cells_apd_4 = np.genfromtxt(cells_apd_filename_4)                       # The cells APD are unordered !!!
    cells_apd_5 = np.genfromtxt(cells_apd_filename_4_1)                       # The cells APD are unordered !!!

    # Get the index and x positions of each cell
    sorted_ids = []
    for i in range(len(sorted_cells_positions)):
        sorted_ids.append(int(sorted_cells_positions[i][0]))
    sorted_x = sorted_cells_positions[:,1]

    mean_apd_column_1 = []
    mean_apd_column_2 = []
    mean_apd_column_3 = []
    mean_apd_column_4 = []
    mean_apd_column_4_1 = []
    x_plot = []
    dx = 100.0
    for i in range(200):
        cells_inside_column = []
    
        minx = i*dx
        maxx = (i+1)*dx

        x_plot.append(i*dx + (dx/2.0))

        print("Getting cells inside interval [%g,%g] ..." % (minx,maxx))
        for j in range(len(sorted_x)):
            if (sorted_x[j] > minx and sorted_x[j] < maxx):
                cells_inside_column.append(sorted_ids[j])
        
        mean_apd_1 = 0.0
        mean_apd_2 = 0.0
        mean_apd_3 = 0.0
        mean_apd_4 = 0.0
        mean_apd_5 = 0.0
        for j in range(len(cells_inside_column)):
            index = cells_inside_column[j]
            apd_1 = cells_apd_1[index]
            apd_2 = cells_apd_2[index]
            apd_3 = cells_apd_3[index]
            apd_4 = cells_apd_4[index]
            apd_5 = cells_apd_5[index]
            #print("%d %g %g %g %g" % (index,sorted_cells_positions[j][1],sorted_cells_positions[j][2],sorted_cells_positions[j][3],cells_apd[index]))

            mean_apd_1 = mean_apd_1 + apd_1
            mean_apd_2 = mean_apd_2 + apd_2
            mean_apd_3 = mean_apd_3 + apd_3
            mean_apd_4 = mean_apd_4 + apd_4
            mean_apd_5 = mean_apd_5 + apd_5
        mean_apd_1 = mean_apd_1 / len(cells_inside_column)
        mean_apd_2 = mean_apd_2 / len(cells_inside_column)
        mean_apd_3 = mean_apd_3 / len(cells_inside_column)
        mean_apd_4 = mean_apd_4 / len(cells_inside_column)
        mean_apd_5 = mean_apd_5 / len(cells_inside_column)
        mean_apd_column_1.append(mean_apd_1)
        mean_apd_column_2.append(mean_apd_2)
        mean_apd_column_3.append(mean_apd_3)
        mean_apd_column_4.append(mean_apd_4)
        mean_apd_column_4_1.append(mean_apd_5)
    
    #print(x_plot)
    #print(mean_apd_column)

    plot_apd_over_a_line(x_plot,mean_apd_column_1,mean_apd_column_2,mean_apd_column_3,mean_apd_column_4,mean_apd_column_4_1)
    

if __name__ == "__main__":
	main()
