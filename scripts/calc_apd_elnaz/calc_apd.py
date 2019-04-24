# Authors: Lucas Berg and Elnaz
# 
# Program that calculates the APD for a single cell in the grid by passing its 
# transmembrane potential over the time as input.
#
# This program also works with multiple AP's !

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
    
    v = vms[start:end]

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

    print("Start = %g -- Finish = %g -- Peak = %d -- Ref = %d" % (start,end,t_peak,t_ref))
    print("APD = %g ms" % (t_ref - t_peak))

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
    if ( len(sys.argv) != 7 ):
        print("-------------------------------------------------------------------------------------------------------------")
        print("Usage:> python %s <input_filename> <num_aps> <period> <ms_each_step>" % sys.argv[0]) 
        print("                  <print_rate> <APD_percentage>")
        print("-------------------------------------------------------------------------------------------------------------")
        print("<input_filename> = Input filename of the simulation")
        print("<num_aps> = Number of AP's to be used for the APD calculation")
        print("<period> = Period of each action potential")
        print("<ms_each_step> = Number of milliseconds from each timestep")
        print("     this value can be calculated as:")
        print("         num_files = (simulation_time) / (dt * print_rate)")
        print("         ms_each_step = (simulation_time) / (num_files)")
        print("     Where the values <simulation_time>, <dt> and <print_rate> are all")
        print("     given in the configuration file of the simulation.")
        print("<print_rate> = The print rate of the simulation")
        print("<APD_percentage> = APD percentage to be used. e.g: APD_90 = 90")
        print("-------------------------------------------------------------------------------------------------------------")
        print("Example:> python calc_apd.py aps/ap-fhn-4aps-cell-1000 4 250 1 100 90")
        print("          python calc_apd.py aps/ap-fhn-cell-1000.txt 1 250 2 100 90")
        print("          python calc_apd.py aps/ap-tt-cell-1000.txt 1 500 2 100 90")
        print("-------------------------------------------------------------------------------------------------------------")
        return 1

    # Get user inputs
    ap_file_name = sys.argv[1]
    num_aps = int(sys.argv[2])
    period = int(sys.argv[3])
    ms_each_step = float(sys.argv[4])
    print_rate = int(sys.argv[5])
    APD_percentage = float(sys.argv[6])

    # Get the transmembrane potential from the input file
    ap_file = open(ap_file_name)
    vms = [float(data) for data in ap_file]
    num_timesteps = len(vms)

    for i in range(num_aps):

        start = i*period
        finish = start + period
    
        print("==================================================================")
        print("AP %d -- APD_percentage = %g" % (i+1,APD_percentage))
        # Call the APD calculation function
        apd = calc_apd(vms,start,finish,ms_each_step,APD_percentage*0.01)
        print("==================================================================")

if __name__ == "__main__":
    main()
