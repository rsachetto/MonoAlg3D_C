# Author: Lucas Berg
# This script only works with one action potential

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

    # Convert from Numpy array to list of float
    d = data[start:(start+end)].tolist()

    max_v = max(d)

    max_index = d.index(max_v) + start

    return max_index


def index_activation(data, start=0):
    d = data[start:]

    for i, v in enumerate(d):
        if d[i + start] < 0.0 < d[i + start + 1]:
            return i+start

def calc_time_reference (vms, start, end, h, refV):
    
    v = vms[start:(start+end)].tolist()

    for i in range(len(v)):
        if (v[i] < refV):
            return (i+start), (i+start)*h
    return -1, -1

def calc_max_min_ref_potential (data,apd_p,start,end):
    
    d = data[start:(start+end)].tolist()

    max_v = max(d)
    min_v = min(d)
    ref_v = min_v + (max_v - min_v)*(1.0 - apd_p)

    return max_v, min_v, ref_v

def calc_max_derivative (vms,start,end,h):

    v = vms[start:(start+end)].tolist()
    n = len(v)
    dvdt = range(0,n)

    for i in range(1,n):
        dvdt[i] = abs(v[i] - v[i-1]/h)
    
    max_dvdt = max(dvdt)

    max_index = dvdt.index(max_dvdt) + start

    return max_index, max_index*h

def calc_ap_limits (vms, start, end):
    
    v = vms[start:(start+end)].tolist()

    maxV = max(vms)
    minV = min(vms)

    index_maxV = v.index(maxV) + start

    return index_maxV, maxV, minV

def calc_reference_potential (minV, maxV, percentage):

    refV = minV + ((maxV - minV)*(1.0 - percentage))

    return refV

def calc_apd (vms, start, end, h, cell_id, percentage=0.9):

    index_maxV, maxV, minV = calc_ap_limits(vms,start,end)
    #print("%d %g %g" % (index_maxV,maxV,minV))

    refV = calc_reference_potential(minV,maxV,percentage)
    #print("%g" % (refV))

    index_peak, t_peak = calc_max_derivative(vms,start,end,h)
    #print("%d %g" % (index_peak,t_peak))

    index_ref, t_ref = calc_time_reference(vms,index_maxV,end,h,refV)
    #print("%d %g" % (index_ref,t_ref))

    if (index_ref == -1):
        print("[Cell %d] ERROR! Could not find reference potential!" % cell_id)
        sys.exit(1)
    
    return (t_ref - t_peak)

def main ():
    if ( len(sys.argv) != 6 ):
        print("-------------------------------------------------------------------------------------------------------------")
        print("Usage:> python %s <output_dir> <num_aps> <ms_each_step> <total_num_cells> <print_rate>" % sys.argv[0])
        print("-------------------------------------------------------------------------------------------------------------")
        print("<output_dir> = Output directory of the simulation")
        print("<total_num_cells> = Total number of cells to calculate the APD")
        print("<num_aps> = Number of AP's to be used for the APD calculation")
        print("<ms_each_step> = Number of milliseconds from each timestep")
        print("     this value can be calculated as:")
        print("         num_files = (simulation_time) / (dt * print_rate)")
        print("         ms_each_step = (simulation_time) / (num_files)")
        print("     Where the values <simulation_time>, <dt> and <print_rate> are all")
        print("     given in the configuration file of the simulation.")
        print("<print_rate> = The print rate of the simulation")
        print("-------------------------------------------------------------------------------------------------------------")
        print("Example:> python calc_full_apd.py ../../outputs/plain_100_100_100_fhn 1 2 10000 100")
        print("-------------------------------------------------------------------------------------------------------------")
        return 1

    # Get user inputs
    output_dir = sys.argv[1]
    num_aps = int(sys.argv[2])
    ms_each_step = float(sys.argv[3])
    total_num_cells = int(sys.argv[4])
    print_rate = int(sys.argv[5])

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
    out_file = open("output/cells-apd.txt","w")

    print("[!] Calculating APD's ...")
    start = time.time()
    # Iterate over each cell
    for i in range(num_cells):
        # Get the transmembrane potential from the current cell
        vms = ap_data[i]

        # Calculates its APD
        apd = calc_apd(vms,0,1000,ms_each_step,i,0.9)

        # Write its value on the output file
        out_file.write("%g\n" % (apd))
    end = time.time()
    print("Elapsed time = %g seconds" % (end - start))
    
    out_file.close()

if __name__ == "__main__":
    main()
