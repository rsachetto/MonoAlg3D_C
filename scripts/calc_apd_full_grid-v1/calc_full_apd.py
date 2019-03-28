import sys
import subprocess
import time

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


# Returns the timestep where the maximum potential occurs
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

def calc_apd (cell_id, num_aps, ms_each_step):
    apds = []
    slope_s = 0
    slope_e = 0
    max_ind = 0

    ap_file = open("aps/cell-%d.txt" % cell_id)

    ap_data = [float(data) for data in ap_file]

    for j in range(num_aps):
            
        #print("-------------------------------------------------------------")
        slope_s = slope_start(ap_data,slope_e,0.5,ms_each_step)
        #print("Slope start = %g ms\n" % (slope_s*ms_each_step))   

        max_ind = max_index(ap_data, slope_e, 1050)
        #print("Max ind = %d\n" % (max_ind))

        slope_e = slope_end(ap_data, max_ind, 0.0005, ms_each_step)
        #print("Slope end = %g\n" % (slope_e*ms_each_step))
        #print("-------------------------------------------------------------")
        
        apds.append((slope_e-slope_s)*ms_each_step)      
        
    apd = sum(apds)/len(apds)
    return apd

def main ():

    if ( len(sys.argv) != 5 ):
        print("-------------------------------------------------------------------------")
        print("Usage:> python %s <output_dir> <total_num_cells> <num_aps> <ms_each_step>" % sys.argv[0])
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
        print("-------------------------------------------------------------------------")
        print("Example:> python calc_full_apd.py ../../outputs/plain_100_100_100_fhn 10 1 2")
        print("-------------------------------------------------------------------------")
        return 1

    output_dir = sys.argv[1]
    total_num_cells = int(sys.argv[2])
    num_aps = int(sys.argv[3])
    ms_each_step = float(sys.argv[4])

    out_file = open("output/cells-apd.txt","w")

    start = time.time()
    # Iterate over each cell
    for i in range(1,20):
    #for i in range(1,total_num_cells):
        print("=============================================================")
        print("Working on cell %d" % i)

        cmd = "./getAps.sh %s V 6 aps/cell-%d.txt %d" % (output_dir,i,i)
        rc = subprocess.call( cmd, shell=True )

        apd = calc_apd(i,num_aps,ms_each_step)

        out_file.write("%g\n" % (apd))
        
        #out_file.write("Cell %d -- APD = %g\n" % (i,apd))
        #print("Cell %d -- APD = %g" % (i,apd))
        print("=============================================================")
    end = time.time()
    print("Elapsed time = %g seconds" % (end - start))
    
    out_file.close()

if __name__ == "__main__":
    main()
