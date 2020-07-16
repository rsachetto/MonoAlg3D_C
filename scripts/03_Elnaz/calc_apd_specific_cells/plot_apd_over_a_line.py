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

def plot_apd_over_a_line (x, apd, f):
	xnew = np.linspace(min(x), max(x), num=1000, endpoint=True)
	plt.grid()
	#plt.plot(x, apd, 'o')
	#plt.plot(x, apd, label="APD", c="black", linewidth=1.0)
	plt.plot(xnew, f(xnew), c="black", linewidth=1.0)
	plt.xlabel("x (um)",fontsize=15)
	plt.ylabel("APD (ms)",fontsize=15)
	plt.title("Action potential duration (APD)",fontsize=14)
	plt.legend(loc=0,fontsize=14)
	#plt.show()
	plt.savefig("ap.pdf")

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

def interpolate_data (x,y):
	#f = interp1d(x, y)
	f = interp1d(x, y, kind='cubic')
	#tck = interpolate.splrep(x, y, s=0)
	#ynew = interpolate.splev(x, tck, der=0)
	
	return f
	#return ynew

def main():
	
	if len(sys.argv) != 3:
		print("-------------------------------------------------------------------------")
		print("Usage:> python %s <cells_positions_filename> <cells_apd_filename>" % sys.argv[0])
		print("-------------------------------------------------------------------------")
		print("<cells_positions_filename> = Input file with the positions of each cell")
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

	xnew = []
	ynew = []
	for i in range(len(x)):
		if (i % 8 == 0):
			xnew.append(x[i])
			ynew.append(cells_apd[i])

	#f = interpolate_data(x,cells_apd)
	f = interpolate_data(xnew,ynew)
	#ynew = interpolate_data(x,cells_apd)

	plot_apd_over_a_line(xnew,ynew,f)
	#plot_apd_over_a_line(x,ynew)

if __name__ == "__main__":
	main()
