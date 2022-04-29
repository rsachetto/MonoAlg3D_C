# ==============================================================================================================================
# Script used to plot the state-vector traces from a simulation that uses the 'save_as_text_or_binary' [save_result].
#  1) Run the simulation and store the 'V_it_<id>.txt' files in the outputs folder
#  2) Define the cell coordinates to get the traces. Save the coordinates in a file
#  3) Write the cell model state-vector names accordindly to the sv[] array from the model
#  4) Run the script to get the traces
# ------------------------------------------------------------------------------------------------------------------------------
# Remarks: 
#   - Multiple cells plot only works for only ventricular simulations
# Author: Lucas Berg
# ==============================================================================================================================
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def get_correct_indexes (sv):
    nlin, ncol = np.shape(sv)

    # Capture the state-vector indexes using Lucas`s scheme 
    v       = sv[:,0]
    CaMKt   = sv[:,1]
    cass    = sv[:,2]
    nai     = sv[:,3]
    nass    = sv[:,4]
    ki      = sv[:,5]
    kss     = sv[:,6]
    cansr   = sv[:,7]
    cajsr   = sv[:,8]
    cai     = sv[:,9]
    m       = sv[:,10]
    h       = sv[:,11]
    j       = sv[:,12]
    hp      = sv[:,13]
    jp      = sv[:,14]
    mL      = sv[:,15]
    hL      = sv[:,16]
    hLp     = sv[:,17]
    a       = sv[:,18]
    iF      = sv[:,19]
    iS      = sv[:,20]
    ap      = sv[:,21]
    iFp     = sv[:,22]
    iSp     = sv[:,23]
    d       = sv[:,24]
    ff      = sv[:,25]
    fs      = sv[:,26]
    fcaf    = sv[:,27]
    fcas    = sv[:,28]
    jca     = sv[:,29]
    ffp     = sv[:,30]
    fcafp   = sv[:,31]
    nca_ss  = sv[:,32]
    nca_i   = sv[:,33]
    C1      = sv[:,34]
    C2      = sv[:,35]
    C3      = sv[:,36]
    I       = sv[:,37]
    O       = sv[:,38]
    xs1     = sv[:,39]
    xs2     = sv[:,40]
    Jrel_np = sv[:,41]
    Jrel_p  = sv[:,42]

    # Change the columns to the Original mapping
    sv_new = np.zeros((nlin,ncol))
    sv_new[:,0]  = v
    sv_new[:,1]  = nai
    sv_new[:,2]  = nass
    sv_new[:,3]  = ki
    sv_new[:,4]  = kss
    sv_new[:,5]  = cai
    sv_new[:,6]  = cass
    sv_new[:,7]  = cansr
    sv_new[:,8]  = cajsr
    sv_new[:,9]  = m
    sv_new[:,10] = hp
    sv_new[:,11] = h
    sv_new[:,12] = j
    sv_new[:,13] = jp
    sv_new[:,14] = mL
    sv_new[:,15] = hL 
    sv_new[:,16] = hLp
    sv_new[:,17] = a
    sv_new[:,18] = iF
    sv_new[:,19] = iS
    sv_new[:,20] = ap
    sv_new[:,21] = iFp
    sv_new[:,22] = iSp
    sv_new[:,23] = d
    sv_new[:,24] = ff
    sv_new[:,25] = fs
    sv_new[:,26] = fcaf
    sv_new[:,27] = fcas
    sv_new[:,28] = jca
    sv_new[:,29] = nca_ss
    sv_new[:,30] = nca_i
    sv_new[:,31] = ffp
    sv_new[:,32] = fcafp
    sv_new[:,33] = xs1
    sv_new[:,34] = xs2
    sv_new[:,35] = Jrel_np
    sv_new[:,36] = CaMKt
    sv_new[:,37] = C1
    sv_new[:,38] = C2
    sv_new[:,39] = C3
    sv_new[:,40] = O
    sv_new[:,41] = I
    sv_new[:,42] = Jrel_p

    return sv_new

def write_limpetgui (t,sv,cell_id):
    nlin, ncol = np.shape(sv)
    filename = "outputs/cell_%d/sv_cell_limpetgui_%d.txt" % (cell_id,cell_id)
    file = open(filename,"w")
    for i in range(nlin):
        file.write("%e " % (t[i]))
        for j in range(ncol):
            file.write("%e " % (sv[i][j]))
        file.write("\n")
    file.close()

def plot_all_state_vector(t, sv_names, num_cells_to_plot):
    for i in range(num_cells_to_plot):
        print("[INFO] Plotting cell %d ..." % (i))
        filename = "outputs/cell_%d/sv_cell_%d.txt" % (i,i)
        sv = np.genfromtxt(filename)
        
        sv = get_correct_indexes(sv)
        write_limpetgui(t,sv,i)

        nlin, ncol = np.shape(sv)
        # Sanity check
        if (ncol != len(sv_names)):
            print("[-] ERROR! Number of columns in the state-vector is different \
                than the state-vector names array")
            print("%d %d" % (ncol,len(sv_names)))
            sys.exit(1)

        for j in range(ncol):
            plt.plot(t,sv[:,j],label="sv[%d] = %s" % (j,sv_names[j]))
            plt.title("%s" % (sv_names[j]))
            plt.legend(loc=0,fontsize=14)
            plt.savefig("outputs/cell_%d/cell_%d__%s.png" % (i,i,sv_names[j]), dpi=150)
            plt.clf()

def read_cells_coordinates (filename):
    cells_coord = []
    file = open(filename,"r")
    for line in file:
        tokens = line.split()
        cells_coord.append([float(tokens[0]),float(tokens[1]),float(tokens[2])])
    file.close()
    return cells_coord

def read_state_vector_file (filename):
    file = open(filename,"r")
    cells_sv = []
    # For each cell
    for line in file:
        cell_sv = []
        tokens = line.split(",")
        # For each state-vector value
        for value in tokens:
            cell_sv.append(float(value))
        cells_sv.append(cell_sv)
    file.close()
    return cells_sv

def read_state_vector_names (filename):
    sv_names = []
    file = open(filename,"r")
    for line in file:
        sv_names.append(line[0:len(line)-1])
    file.close()
    return sv_names

def write_state_vector (cells_to_plot,cells_sv):
    counter = 0
    for cell_to_plot in cells_to_plot:
        filename = "outputs/cell_%d/sv_cell_%d.txt" % (counter,counter)
        file = open(filename,"a")   # Append the results
        for cell_sv in cells_sv:
            if (cell_to_plot[0] == cell_sv[0] and cell_to_plot[1] == cell_sv[1] and cell_to_plot[2] == cell_sv[2]):
                for i in range(6,len(cell_sv)):
                    file.write("%g " % (cell_sv[i]))
                file.write("\n")
        counter = counter + 1
        file.close()

def create_directory_and_check_old_files (num_cells_to_plot):
    try:
        print("[INFO] Removing old files ...")
        os.system("rm -r outputs/cell_*")
    except OSError as error: 
        print(error)
    
    try: 
        os.mkdir("outputs") 
    except OSError as error: 
        print(error)
    for i in range(num_cells_to_plot):
        try: 
            os.mkdir("outputs/cell_%d" % (i)) 
        except OSError as error: 
            print(error)

def main():
	
    if len(sys.argv) != 7:
        print("----------------------------------------------------------------------------------------------------------------------")
        print("Usage:> python %s <input_folder> <tmax> <dt> <printrate> <cells_coordinates_file> <sv_name_array>" % sys.argv[0])
        print("----------------------------------------------------------------------------------------------------------------------")
        print("<input_file> = Input folder with the state-vector data")
        print("<tmax> = Maximum simulation time")
        print("<dt> = Timestep used in the simulation")
        print("<printrate> = Print rate used in the simulation")
        print("<cells_coordinates_file> = File with the cell coordinates to print")
        print("<sv_name_array> = State-vector variable name array")
        print("----------------------------------------------------------------------------------------------------------------------")
        return 1

    input_folder = sys.argv[1]
    tmax = float(sys.argv[2])
    dt = float(sys.argv[3])
    printrate = int(sys.argv[4])
    cell_coordinates_filename = sys.argv[5]
    sv_names_filename = sys.argv[6]

    # Get the cells coordinates to plot
    cells_to_plot = read_cells_coordinates(cell_coordinates_filename)
    num_cells_to_plot = len(cells_to_plot)

    # Create output folder
    create_directory_and_check_old_files(num_cells_to_plot)

    # Get state-vector names
    sv_names = read_state_vector_names(sv_names_filename)
    
    # Total number of generated files by the simulation
    num_files = int(tmax / (dt*printrate))
    t = [i*dt*printrate for i in range(num_files)]
    t = np.linspace(0,1000,50000)

    #for i in range(num_files):
    for i in range(4950000,5000000):
        # Calculate the current iteration number
        step_id = i*printrate

        # Open the current iteration results file
        filename = "%s/V_it_%d.txt" % (input_folder,step_id)

        # Get the state-vector for the current iteration
        cells_sv = read_state_vector_file(filename)
        neq = len(cells_sv)-7

        # Write the state-vector for each plot cell by appending the lines
        write_state_vector(cells_to_plot,cells_sv)

    # Call the plot function
    plot_all_state_vector(t,sv_names,num_cells_to_plot)

if __name__ == "__main__":
	main()
