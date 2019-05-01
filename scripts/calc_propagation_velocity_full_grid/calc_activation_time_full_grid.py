# Author: Lucas Berg
#
# Program that calculates the APD of every cell in the grid by passing a certain percentage. 
#   e.g: APD_90 --> percentage = 90
#        APD_80 --> percentage = 80
#        APD_70 --> percentage = 70
#
# This program also works with multiples AP's !

import sys
import subprocess
import time
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

def calc_ap_limits (vms, start, end):
    
    v = vms[start:end].tolist()

    maxV = max(v)
    minV = min(v)

    index_maxV = v.index(maxV) + start

    return index_maxV, maxV, minV

def calc_reference_potential (minV, maxV, percentage):

    refV = minV + (maxV - minV)*(1.0 - percentage)

    return refV        

def calc_time_reference (vms, start, end, h, refV):
    
    v = vms[start:(start+end)]

    for i in range(len(v)):
        if (v[i] < refV):
            return (i+start), (i+start)*h
    return -1, -1

def calc_max_min_ref_potential (data,apd_p,start,end):
    
    d = data[start:end]

    max_v = max(d)
    min_v = min(d)
    ref_v = min_v + (max_v - min_v)*(1.0 - apd_p)

    return max_v, min_v, ref_v

def calc_max_derivative (vms,start,end,h):
    v = vms[start:end]
    n = len(v)
    dvdt = range(0,n)

    for i in range(1,n):
        dvdt[i] = abs(v[i] - v[i-1]/h)
    
    max_dvdt = max(dvdt)

    max_index = dvdt.index(max_dvdt) + start

    return max_index, max_index*h

def calc_activation_time (vms, h, start, end):

    index_peak, t_peak = calc_max_derivative(vms,start,end,h)
    
    return (t_peak)

def main ():
    if ( len(sys.argv) != 5 ):
        print("-------------------------------------------------------------------------------------------------------------")
        print("Usage:> python %s <output_dir> <ms_each_step> <total_num_cells> <print_rate>" % sys.argv[0])
        print("-------------------------------------------------------------------------------------------------------------")
        print("<output_dir> = Output directory of the simulation")
        print("<ms_each_step> = Number of milliseconds from each timestep")
        print("     this value can be calculated as:")
        print("         num_files = (simulation_time) / (dt * print_rate)")
        print("         ms_each_step = (simulation_time) / (num_files)")
        print("     Where the values <simulation_time>, <dt> and <print_rate> are all")
        print("     given in the configuration file of the simulation.")
        print("<total_num_cells> = Total number of cells to calculate the APD")
        print("<print_rate> = The print rate of the simulation")
        print("-------------------------------------------------------------------------------------------------------------")
        print("Example:> python calc_activation_time_full_grid.py ../../outputs/plain_100_100_100_fhn 2 10000 100")
        print("          python calc_activation_time_full_grid.py ../../outputs/plain_100_100_100_tentusscher 2 10000 100")
        print("-------------------------------------------------------------------------------------------------------------")
        return 1

    # Get user inputs
    output_dir = sys.argv[1]
    ms_each_step = float(sys.argv[2])
    total_num_cells = int(sys.argv[3])
    print_rate = int(sys.argv[4])

    # Get the transmembrane potential for all timesteps and every single cell in the grid
    print("[!] Calling 'getAps.sh' script ...")
    cmd = "./getAps.sh %s V 6 %d %d" % (output_dir,total_num_cells,print_rate)
    rc = subprocess.call( cmd, shell=True )

    # Open the generated file from the previous script as a Numpy array and get its dimensions
    timesteps_file = open("timesteps.txt")
    ap_data = np.genfromtxt(timesteps_file)
    num_cells = ap_data.shape[0]
    num_timesteps = ap_data.shape[1]
    
    # Open the output file
    out_file = open("inputs/cells-at.txt","w")

    print("[!] Calculating Activation Times ...")
    
    begin = time.time()
    # Iterate over each cell
    for i in range(num_cells):
        # Get the transmembrane potential from the current cell
        vms = ap_data[i]

        # Calculate the activation time of that cell
        at = calc_activation_time(vms,ms_each_step,0,500)

        # Write the value to an output file
        out_file.write("%g\n" % (at))
    end = time.time()
    
    print("[!] Elapsed time = %g seconds" % (end - begin))
    
    out_file.close()
    print("[+] Output file saved in 'inputs/cell-at.txt'")

if __name__ == "__main__":
    main()
