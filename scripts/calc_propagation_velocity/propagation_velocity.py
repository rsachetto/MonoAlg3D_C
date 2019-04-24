# Author: Lucas Berg
# ---------------------------------------------------------------------------------------------------
# GUIDE:
#   To generate the input files we need to use the "getAps.sh" and "cell_positions_calculator" programs.
#   The "getAps.sh" will compute the transmembrane potential from the two reference cells set for
#  the propagation velocity calculus. This way we could compute their activation times using the
#  "slope_start" function.
#   The "cell_positions_calculator" is needed because the coordinates from the two reference cells
#  must be known, in order to calculate the distance between them.
# ---------------------------------------------------------------------------------------------------

import sys
import numpy as np

def forwarddiff(y, h):
    n = len(y)
    res = []

    for i in range(1, n):
        res.append((y[i] - y[i-1]) / h)

    return res


def slope_start(data, start=0, epsilon=0.0001, h=1.0):

    d = data[start:]
    n = len(d)

    for i in range(1,n):
        if abs(d[i] - d[i-1]/h) > epsilon:
            return i+start


def slope_end(data, start=0, epsilon=0.0001, h=1.0):

    d = data[start:]
    n = len(d)

    for i in range(1,n):
        if abs(d[i] - d[i-1]/h) < epsilon:
            return i+start

def max_index(data, start, end):

    d = data[start:(start+end)]

    max_v = max(d)

    max_index = d.index(max_v) + start

    return max_index


def index_activation(data, start=0):
    d = data[start:]

    for i, v in enumerate(d):
        if d[i + start] < 0.0 < d[i + start + 1]:
            return i+start

def read_ap_data (ap_filename):

    ap_file = open(ap_filename)

    ap_data = [float(data) for data in ap_file]

    return ap_data

def read_cells_positions (cells_filename):
    
    cells_data = np.genfromtxt(cells_filename)

    return cells_data

def calculate_cell_distance (cells_positions, cell_1_id, cell_2_id):
    
    x1, y1, z1 = cells_positions[cell_1_id]
    x2, y2, z2 = cells_positions[cell_2_id]

    d = np.sqrt( np.power(x2-x1,2) + np.power(y2-y1,2) + np.power(z2-z1,2) )

    return d

def calculate_velocity (d,t1,t2):
    
    v = d / abs(t2-t1)

    return v

def main ():

    if ( len(sys.argv) != 6 ):
        print("-------------------------------------------------------------------------------------------------------")
        print("Usage:> python %s <ap_file_name_cell_1> <ap_file_name_cell_2> <cell_positions_file_name>" % sys.argv[0])
        print("                  <cell_1_id> <cell_2_id>")
        print("-------------------------------------------------------------------------------------------------------")
        print("<ap_file_name_cell_1> = Input file with the AP's from cell 1 at each timestep")
        print("<ap_file_name_cell_2> = Input file with the AP's from cell 2 at each timestep")
        print("<cell_positions_file_name> = Input file with the cell positions")
        print("<cell_1_id> = Index of the first reference cell")
        print("<cell_2_id> = Index of the second reference cell")
        print("-------------------------------------------------------------------------------------------------------")
        return 1
    
    ap_file_name_cell_1 = sys.argv[1]
    ap_file_name_cell_2 = sys.argv[2]
    cell_positions_file_name = sys.argv[3]
    cell_1_id = int(sys.argv[4])
    cell_2_id = int(sys.argv[5])

    ap_data_cell_1 = read_ap_data(ap_file_name_cell_1)
    ap_data_cell_2 = read_ap_data(ap_file_name_cell_2)
    cells_positions = read_cells_positions(cell_positions_file_name)
    
    t1 = slope_start(ap_data_cell_1)
    t2 = slope_start(ap_data_cell_2)

    d = calculate_cell_distance(cells_positions,cell_1_id,cell_2_id)
    
    v = calculate_velocity(d,t1,t2)

    print("Activation time of cell 1 = %g ms" % t1)
    print("Activation time of cell 2 = %g ms" % t2)
    print("Distance = %g cm" % d)
    print("Velocity = %g cm/ms" % v)

if __name__ == "__main__":
    main()
