# ============================================================================================================
# Author: Lucas Berg
#
# Program that calculates the APD of specific cells in the grid by passing a file with indexes
# from the desired cells and a certain percentage for the APD. 
#   e.g: APD_90 --> percentage = 90
#        APD_80 --> percentage = 80
#        APD_70 --> percentage = 70
#
# This program also works with multiples AP's !
#
# IMPORTANT: This program must be used together 'cell_position_calculator_over_a_line' and only works
#           with a plain mesh. 
# ============================================================================================================

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

def calc_apd (vms, start, end, h, percentage):

    index_maxV, maxV, minV = calc_ap_limits(vms,start,end)

    refV = calc_reference_potential(minV,maxV,percentage)

    index_peak, t_peak = calc_max_derivative(vms,start,end,h)

    index_ref, t_ref = calc_time_reference(vms,index_maxV,end,h,refV)

    #print("Start = %g -- Finish = %g -- Peak = %d -- Ref = %d" % (start,end,t_peak,t_ref))
    #print("APD = %g ms" % (t_ref - t_peak))

    if (index_ref == -1):
        print("[-] ERROR! Could not find reference potential!")
        sys.exit(1)
    
    return (t_ref - t_peak)
    

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

def main ():
    if ( len(sys.argv) != 9 ):
        print("-------------------------------------------------------------------------------------------------------------")
        print("Usage:> python %s <output_dir> <cell_positions_filename>" % sys.argv[0])
        print("                  <num_aps> <period> <ms_each_step> <total_num_cells>")
        print("                  <print_rate> <APD_percentage>")
        print("-------------------------------------------------------------------------------------------------------------")
        print("<output_dir> = Output directory of the simulation")
        print("<cell_positions_filename> = Filename with the indexes and positions from the specific cells")
        print("<num_aps> = Number of AP's to be used for the APD calculation")
        print("<period> = Period of each action potential occurs")
        print("<ms_each_step> = Number of milliseconds from each timestep")
        print("     this value can be calculated as:")
        print("         num_files = (simulation_time) / (dt * print_rate)")
        print("         ms_each_step = (simulation_time) / (num_files)")
        print("     Where the values <simulation_time>, <dt> and <print_rate> are all")
        print("     given in the configuration file of the simulation.")
        print("<total_num_cells> = Total number of cells to calculate the APD")
        print("<print_rate> = The print rate of the simulation")
        print("<APD_percentage> = Percentage to be used as reference on the APD calculus")
        print("-------------------------------------------------------------------------------------------------------------")
        print("Example:> python calc_apd_specific_cells.py ../../outputs/plain_mixed_models_tt inputs/cells_positions_inside_region.txt 1 500 1 10000 50 90")
        print("-------------------------------------------------------------------------------------------------------------")
        return 1

    # Get user inputs
    output_dir = sys.argv[1]
    cell_positions_filename = sys.argv[2]
    num_aps = int(sys.argv[3])
    period = int(sys.argv[4])
    ms_each_step = float(sys.argv[5])
    total_num_cells = int(sys.argv[6])
    print_rate = int(sys.argv[7])
    APD_percentage = float(sys.argv[8])

    # Get the transmembrane potential for all timesteps and every single cell in the grid
    print("[!] Calling 'getAps.sh' script ...")
    cmd = "./getAps.sh %s V 6 %d %d" % (output_dir,total_num_cells,print_rate)
    rc = subprocess.call( cmd, shell=True )

    # Open the generated file from the previous script as a Numpy array and get its dimensions
    timesteps_file = open("timesteps.txt")
    ap_data = np.genfromtxt(timesteps_file)
    num_cells = ap_data.shape[0]
    num_timesteps = ap_data.shape[1]

    # Open the input cell positions file
    cell_positions_file = open(cell_positions_filename)
    index_positions = np.genfromtxt(cell_positions_file)
    num_indexes = index_positions.shape[0]
    num_cols = index_positions.shape[1]

    # Open the output file
    out_file = open("output/cells-apd-inside-region.txt","w")

    print("[!] Calculating APD's ...")
    
    begin = time.time()
    # Iterate over each cell inside the region
    for i in range(num_indexes):
        
        # Get the index of the cell
        index = int(index_positions[i][0])

        # Get the transmembrane potential from the current cell
        vms = ap_data[index]

        apds = []
        for j in range(num_aps):
            
            start = j*period
            finish = start + period

            # Calculates its APD
            apd = calc_apd(vms,start,finish,ms_each_step,APD_percentage*0.01)
            
            apds.append(apd)

        # Write the mean APD value on the output file
        out_file.write("%g\n" % (sum(apds)/len(apds)))
    end = time.time()
    
    print("[!] Elapsed time = %g seconds" % (end - begin))
    
    out_file.close()
    print("[+] Output file saved in 'output/cells-apd-inside-region.txt'")

if __name__ == "__main__":
    main()
