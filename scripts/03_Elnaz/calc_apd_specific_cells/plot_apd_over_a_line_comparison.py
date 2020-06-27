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

def plot_apd_over_a_line (x,y1,y2,y3,y4,y4_1):
	
    rms_12 = calculate_rms(y1,y2)
    rms_13 = calculate_rms(y1,y3)
    rms_14 = calculate_rms(y1,y4)
    rms_14_1 = calculate_rms(y1,y4_1)

    plt.grid()
    plt.scatter(x, y1, label="SC1", s=5)
    plt.scatter(x, y2, label="SC2", s=5)
    plt.xlabel("x (um)",fontsize=15)
    plt.ylabel("APD (ms)",fontsize=15)
    plt.title("RMSPE = %g" % (rms_12),fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.savefig("apd-comparison-sc1-sc2.png")
    plt.clf()

    plt.grid()
    plt.scatter(x, y1, label="SC1", s=5)
    plt.scatter(x, y3, label="SC3", s=5)
    plt.xlabel("x (um)",fontsize=15)
    plt.ylabel("APD (ms)",fontsize=15)
    plt.title("RMSPE = %g" % (rms_13),fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.savefig("apd-comparison-sc1-sc3.png")
    plt.clf()

    plt.grid()
    plt.scatter(x, y1, label="SC1", s=5)
    plt.scatter(x, y4, label="SC4", s=5)
    plt.xlabel("x (um)",fontsize=15)
    plt.ylabel("APD (ms)",fontsize=15)
    plt.title("RMSPE = %g" % (rms_14),fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.savefig("apd-comparison-sc1-sc4.png")
    plt.clf()

    plt.grid()
    plt.scatter(x, y1, label="SC1", s=5)
    plt.scatter(x, y4_1, label="SC4_1", s=5)
    plt.xlabel("x (um)",fontsize=15)
    plt.ylabel("APD (ms)",fontsize=15)
    plt.title("RMSPE = %g" % (rms_14_1),fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.savefig("apd-comparison-sc1-sc4_1.png")
    plt.clf()

def sort_apd_by_x (x,apd):
	n = x.shape[0]

	# TODO: Change this bubblesort ...
	for i in range(n):
		for j in range(n):
			if (x[i] < x[j]):
				aux = x[i]
				x[i] = x[j]
				x[j] = aux

				aux = apd[i]
				apd[i] = apd[j]
				apd[j] = aux

def get_data (cells_filename,apd_filename):

    um_to_cm = 1.0e-04

    # Read the input files as Numpy arrays 
    cells_indexes_positions = np.genfromtxt(cells_filename)
    cells_apd = np.genfromtxt(apd_filename)

    # Get the x position from each cell
    x = cells_indexes_positions[:,1]*um_to_cm

    # TODO: Find a way to store (x,cells_apd) as a structure and sort it by the 'x' value
    sort_apd_by_x(x,cells_apd)

    return x, cells_apd

def main():
	
    if len(sys.argv) != 7:
        print("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
        print("Usage:> python %s <cells_positions_filename> <cells_apd_filename_0> <cells_apd_filename_1> <cells_apd_filename_2> <cells_apd_filename_3> <cells_apd_filename_3_1>" % sys.argv[0])
        print("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
        print("<cells_positions_filename> = Input file with the positions of each cell")
        print("<cells_apd_filename_0> = Input file with the APD of each cell")
        print("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
        print("Example:> python %s inputs/cells_positions_inside_region.txt outputs/cells-apd-inside-region.txt" % sys.argv[0])
        return 1

    cells_position_filename = sys.argv[1]
    cells_apd_filename_1 = sys.argv[2]
    cells_apd_filename_2 = sys.argv[3]
    cells_apd_filename_3 = sys.argv[4]
    cells_apd_filename_4 = sys.argv[5]
    cells_apd_filename_4_1 = sys.argv[6]

    x, y1 = get_data(cells_position_filename,cells_apd_filename_1)
    x, y2 = get_data(cells_position_filename,cells_apd_filename_2)
    x, y3 = get_data(cells_position_filename,cells_apd_filename_3)
    x, y4 = get_data(cells_position_filename,cells_apd_filename_4)
    x, y4_1 = get_data(cells_position_filename,cells_apd_filename_4_1)

    plot_apd_over_a_line(x,y1,y2,y3,y4,y4_1)

if __name__ == "__main__":
	main()
