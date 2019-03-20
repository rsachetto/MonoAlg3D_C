import sys

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

def main ():

    if ( len(sys.argv) != 4 ):
        print("-------------------------------------------------------------------------")
        print("Usage:> python %s <ap_file_name> <num_aps> <ms_each_step>" % sys.argv[0])
        print("-------------------------------------------------------------------------")
        print("<ap_file_name> = Input file with the AP's from each timestep")
        print("<num_aps> = Number of AP's to be used for the APD calculation")
        print("<ms_each_step> = Number of milliseconds from each timestep")
        print("     this value can be calculated as:")
        print("         num_files = (simulation_time) / (dt * print_rate)")
        print("         ms_each_step = (simulation_time) / (num_files)")
        print("     Where the values <simulation_time>, <dt> and <print_rate> are all")
        print("     given in the configuration file of the simulation.")
        print("-------------------------------------------------------------------------")
        return 1

    ap_file_name = sys.argv[1]
    num_aps = int(sys.argv[2])
    ms_each_step = float(sys.argv[3])

    ap_file = open(ap_file_name)

    rests = []

    epsilon = 0.003

    ap_data = [float(data) for data in ap_file]

    i = 0

    while True:

        try:
            data1 = ap_data[i]
            data2 = ap_data[i+1]
            i = i + 1

            if abs(data1-data2) < epsilon:
                rests.append(data1)
                rests.append(data2)

        except:
            break

    rest = sum(rests)/len(rests)

    print("Rest: ", rest)
    epsilon = 0.05
    apds = []

    slope_e = 0
    max_ind = 0

    # These settings are working for the Tentusscher -> "elnaz_plain_mesh_tentusscher.ini"
    for j in range(num_aps):

        slope_s = slope_start(ap_data, slope_e, 50)
        print("Slope start = %g\n" % slope_s)

        max_ind = max_index(ap_data, slope_e, 1050)
        print("Max ind = %d\n" % max_ind)

        slope_e = slope_end(ap_data, max_ind+3, 0.005, ms_each_step)
        print("Slope end = %g\n" % slope_e)

        apds.append((slope_e-slope_s)*ms_each_step)

    print("APD: ", sum(apds)/len(apds))

if __name__ == "__main__":
    main()
