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

def sort_apd_by_x (ids,x,y,z):
    # Get the x coordinate from all the cells
    n = x.shape[0]

    for i in range(n):
        for j in range(n):
            if (x[i] < x[j]):
                aux = ids[i]
                ids[i] = ids[j]
                ids[j] = aux

                aux = x[i]
                x[i] = x[j]
                x[j] = aux

                aux = y[i]
                y[i] = y[j]
                y[j] = aux

                aux = z[i]
                z[i] = z[j]
                z[j] = aux
    
def main():
	
    if len(sys.argv) != 2:
        print("-------------------------------------------------------------------------")
        print("Usage:> python %s <cells_positions_filename>" % sys.argv[0])
        print("-------------------------------------------------------------------------")
        print("<cells_positions_filename> = Input file with the positions of each cell")
        print("-------------------------------------------------------------------------")
        print("Example:> python %s inputs/cells_positions_inside_region.txt" % sys.argv[0])
        return 1

    cells_position_filename = sys.argv[1]

    # Read the input files as Numpy arrays
    #print("Reading inputs files ...") 
    cells = np.genfromtxt(cells_position_filename)

    ids = []
    for i in range(len(cells)):
        ids.append(int(cells[i][0]))

    x = cells[:,1]
    y = cells[:,2]
    z = cells[:,3]

    # Sort all the cells using the x position
    #print("Sorting the cells by the x position ...")
    sort_apd_by_x(ids,x,y,z)

    for i in range(len(ids)):
        # Output format: [id x y z APD]
        print("%d %g %g %g" % (ids[i],x[i],y[i],z[i]))
    

if __name__ == "__main__":
	main()
