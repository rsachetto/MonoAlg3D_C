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

def calc_apd (vms, num_aps, ms_each_step):
    apds = []
    slope_s = 0
    slope_e = 0
    max_ind = 0

    for j in range(num_aps):
            
        #print("-------------------------------------------------------------")
        slope_s = slope_start(vms,slope_e,0.5,ms_each_step)
        #print("Slope start = %g ms\n" % (slope_s*ms_each_step))   

        max_ind = max_index(vms, slope_e, 1050)
        #print("Max ind = %d\n" % (max_ind))

        slope_e = slope_end(vms, max_ind, 0.0005, ms_each_step)
        #print("Slope end = %g\n" % (slope_e*ms_each_step))
        #print("-------------------------------------------------------------")
        
        apds.append((slope_e-slope_s)*ms_each_step)      
        
    apd = sum(apds)/len(apds)
    return apd

def clean_files ():
    

def main ():

    if ( len(sys.argv) != 5 ):
        print("-------------------------------------------------------------------------")
        print("Usage:> python %s <output_dir> <num_aps> <ms_each_step> <print_rate>" % sys.argv[0])
        print("-------------------------------------------------------------------------")
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
        print("-------------------------------------------------------------------------")
        print("Example:> python calc_full_apd.py ../../outputs/plain_100_100_100_fhn 1 2 100")
        print("-------------------------------------------------------------------------")
        return 1

    # Get user inputs
    output_dir = sys.argv[1]
    num_aps = int(sys.argv[2])
    ms_each_step = float(sys.argv[3])
    print_rate = int(sys.argv[4])

    # Get the transmembrane potential for all timesteps and every single cell in the grid
    print("[!] Calling 'getAps.sh' script ...")
    cmd = "./getAps.sh %s V 6 %d" % (output_dir,print_rate)
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
        apd = calc_apd(vms,num_aps,ms_each_step)

        # Write its value on the output file
        out_file.write("%g\n" % (apd))
    end = time.time()
    print("Elapsed time = %g seconds" % (end - start))
    
    out_file.close()

if __name__ == "__main__":
    main()
