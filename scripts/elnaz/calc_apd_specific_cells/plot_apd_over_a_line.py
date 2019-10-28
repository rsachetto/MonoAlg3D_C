# =================================================================
# Author: Lucas Berg
#
# Program that plot the APD over a line
# This script only works if the grid is a plain tissue
# =================================================================

import sys
import numpy as np
import matplotlib.pyplot as plt

def plot_apd_over_a_line (x, apd):
	plt.grid()
	plt.plot(x, apd, label="APD", c="black", linewidth=1.0)
	plt.xlabel("x (um)",fontsize=15)
	plt.ylabel("APD (ms)",fontsize=15)
	plt.title("Action potential duration (APD)",fontsize=14)
	plt.legend(loc=0,fontsize=14)
	plt.show()
	#plt.savefig("ap.pdf")

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

def main():
	
	if len(sys.argv) != 3:
		print("-------------------------------------------------------------------------")
		print("Usage:> python %s <cells_positions_filename> <cells_apd_filename>" % sys.argv[0])
		print("-------------------------------------------------------------------------")
		print("<cells_positions_filename> = Input file with the position of each cell")
		print("<cells_apd_filename> = Input file with the APD of each cell")
		print("-------------------------------------------------------------------------")
		print("Example:> python %s inputs/cells_positions_inside_region.txt outputs/cells-apd-inside-region.txt" % sys.argv[0])
		return 1

	cells_position_filename = sys.argv[1]
	cells_apd_filename = sys.argv[2]

	# Read the input files as Numpy arrays 
	cells_indexes_positions = np.genfromtxt(cells_position_filename)
	cells_apd = np.genfromtxt(cells_apd_filename)

	# Get the x position from each cell
	x = cells_indexes_positions[:,1]*1.0e-04

	# TODO: Find a way to store (x,cells_apd) as a structure and sort it by the 'x' value
	sort_apd_by_x(x,cells_apd)

	plot_apd_over_a_line(x,cells_apd)

if __name__ == "__main__":
	main()
