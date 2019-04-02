# Author: Lucas Berg
#
#   Butterfly plotter:
#
# Program that plots the transmembrane potential of every cell in the grid

import sys
import subprocess
import time
import numpy as np
import matplotlib.pyplot as plt

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

def calc_timesteps (num_files,ms_each_step):
    
    t = np.arange(0,num_files)*ms_each_step
    return t

def main ():
    if ( len(sys.argv) != 6 ):
        print("-------------------------------------------------------------------------------------------------------------")
        print("Usage:> python %s <output_dir_1> <output_dir_2> <ms_each_step> <print_rate> <total_num_cells>" % sys.argv[0])
        print("-------------------------------------------------------------------------------------------------------------")
        print("<output_dir_1> = Output directory of the simulation number 1")
        print("<output_dir_2> = Output directory of the simulation number 2")
        print("<ms_each_step> = Number of milliseconds from each timestep")
        print("     this value can be calculated as:")
        print("         num_files = (simulation_time) / (dt * print_rate)")
        print("         ms_each_step = (simulation_time) / (num_files)")
        print("     Where the values <simulation_time>, <dt> and <print_rate> are all")
        print("     given in the configuration file of the simulation.")
        print("<total_num_cells> = Total number of cells to calculate the APD")
        print("<print_rate> = The print rate of the simulation")
        print("-------------------------------------------------------------------------------------------------------------")
        print("Example:> python butterfly_plot.py ../../outputs/plain_100_100_100_fhn 2 100 10000")
        print("-------------------------------------------------------------------------------------------------------------")
        return 1

    # Get user inputs
    output_dir_1 = sys.argv[1]
    output_dir_2 = sys.argv[2]
    ms_each_step = float(sys.argv[3])
    print_rate = int(sys.argv[4])
    total_num_cells = int(sys.argv[5])

    # Get the transmembrane potential for all timesteps and every single cell in the grid
    print("[!] Calling 'getAps.sh' script ...")
    
    print("\t[!] Getting simulation 1 AP's ...")
    cmd = "./getAps.sh %s V 6 %d %d 1" % (output_dir_1,total_num_cells,print_rate)
    rc = subprocess.call( cmd, shell=True )

    print("\t[!] Getting simulation 2 AP's ...")
    cmd = "./getAps.sh %s V 6 %d %d 2" % (output_dir_2,total_num_cells,print_rate)
    rc = subprocess.call( cmd, shell=True )

    # Open the generated files from the previous script as a Numpy array and get its dimensions
    timesteps_file_1 = open("timesteps-1.txt")
    ap_data_1 = np.genfromtxt(timesteps_file_1)

    timesteps_file_2 = open("timesteps-2.txt")
    ap_data_2 = np.genfromtxt(timesteps_file_2)
    
    num_cells = ap_data_1.shape[0]
    num_files = ap_data_1.shape[1]    
    t = calc_timesteps(num_files,ms_each_step)

    print("[!] Plotting AP's ...")
    
    plt.grid()
    plt.xlabel("t (ms)",fontsize=15)
    plt.ylabel("V (mV)",fontsize=15)
    plt.title("Action potential",fontsize=14)
    plt.legend(loc=2,fontsize=14)

    # Iterate over the 2 simulations
    for k in range(2):
        
        # Iterate over each cell
        for i in range(0,10000,100):
            
            # Get the transmembrane potential from the current simulation for the current cell
            if (k == 0):
                vms = ap_data_1[i]
            else:
                vms = ap_data_2[i]
            
            # Plot the AP of the cell using the simulation number to distinguish its color
            if (k == 0):
                plt.plot(t, vms, label="Vm", c="blue", linewidth=1.0)
            else:
                plt.plot(t, vms, label="Vm", c="red", linewidth=1.0)
	
    #plt.show()
    plt.savefig("output/butterfly.pdf")
    print("[+] Output file saved at 'output/butterfly.pdf'")


if __name__ == "__main__":
    main()
