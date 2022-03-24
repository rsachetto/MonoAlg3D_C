# ==============================================================================================================================
# Script used to plot comparison the state-vector traces from two distinct simulations that uses the 'save_as_text_or_binary' 
# [save_result].
#  1) Run the simulations and store the 'V_it_<id>.txt' files in the outputs folder
#  2) Define the cell coordinates to get the traces. Save the coordinates in a file
#  3) Write the cell model state-vector names accordindly to the sv[] array from the model
#  4) Run the script to get the traces
# ------------------------------------------------------------------------------------------------------------------------------
# Remarks: 
#   - Multiple cells plot only works for only ventricular simulations
# Author: Lucas Berg
# ==============================================================================================================================

import sys
import numpy as np
import matplotlib.pyplot as plt

def read_state_variable(filename):
    print("Reading state-vector from '%s' ..." % (filename))
    sv = np.genfromtxt(filename)
    return sv

def plot_state_variable_comparison(t, sv1, sv2, sv_names):

    nlin1, ncol1 = np.shape(sv1)
    nlin2, ncol2 = np.shape(sv2)
    if (ncol1 != len(sv_names) or ncol2 != len(sv_names)):
        print("[-] ERROR! Number of columns in the state-vector is different \
                than the state-vector names array")
        sys.exit(1)

    for j in range(ncol1):
            print("Working on state-variable '%s' ..." % (sv_names[j]))
            #plt.plot(t,sv1[:,j],label="EulerAdapt")
            #plt.plot(t,sv2[:,j],label="RLAdapt")
            #plt.plot(t[245000:],sv1[245000:,j],label="reltol=1e-3")
            #plt.plot(t[245000:],sv2[245000:,j],label="reltol=1e-12")
            plt.plot(t[245000:],sv1[245000:,j],label="dt=EulerAdapt",linewidth=1.0)
            plt.plot(t[245000:],sv2[245000:,j],label="dt=RLAdapt",linewidth=1.0)
            plt.title("SV variable - %s" % (sv_names[j]))
            plt.legend(loc=0,fontsize=14)
            plt.savefig("outputs/comparison_cell__%s.png" % (sv_names[j]), dpi=150)
            plt.clf()

def read_state_vector_names (filename):
    sv_names = []
    file = open(filename,"r")
    for line in file:
        sv_names.append(line[0:len(line)-1])
    file.close()
    return sv_names

def calc_error (sv_ref, sv_aprox, sv_names):
    nlin_ref, ncol_ref = np.shape(sv_ref)
    nlin_aprox, ncol_aprox = np.shape(sv_aprox)
    # Sanity check
    if (nlin_ref != nlin_aprox or ncol_ref != ncol_aprox):
        print("[-] ERROR! The number of lines or coluns between the state-vectors are different!")
        sys.exit(1)

    sv_ref = sv_ref[245000:,:]
    sv_aprox = sv_aprox[245000:,:]
    rel_error = np.abs(sv_ref-sv_aprox)/np.abs(sv_ref+1e-16)
    rel_error = rel_error.sum(axis=0) / nlin_ref * 100.0

    for i in range(ncol_ref):
        print("sv[%d]{%s} = %.2lf" % (i,sv_names[i],rel_error[i]))

    return rel_error

def main():
    if len(sys.argv) != 7:
        print("----------------------------------------------------------------------------------------------------------------------------------")
        print("Usage:> python %s <input_file_1> <input_file_2> <tmax> <dt> <print_rate> <sv_name>" % sys.argv[0])
        print("----------------------------------------------------------------------------------------------------------------------------------")
        print("<input_file_1> = Input file with the state-vector from the first simulation")
        print("<input_file_2> = Input file with the state-vector from the second simulation")
        print("<tmax> = Timestep value used for the simulation")
        print("<dt> = Timestep value used for the simulation")
        print("<print_rate> = Print rate used for the simulation")
        print("<sv_name_array> = State-vector variable name array")
        print("----------------------------------------------------------------------------------------------------------------------------------")
        return 1

    input_file_1 = sys.argv[1]
    input_file_2 = sys.argv[2]
    tmax = float(sys.argv[3])
    dt = float(sys.argv[4])
    printrate = int(sys.argv[5])
    sv_names_filename = sys.argv[6]

    # Total number of generated files by the simulation
    num_files = int(tmax / (dt*printrate)+1)
    t = [i*dt*printrate for i in range(num_files)]

    # Read state-vector values
    sv1 = read_state_variable(input_file_1)
    sv2 = read_state_variable(input_file_2)

    # Get state-vector names
    sv_names = read_state_vector_names(sv_names_filename)

    #calc_error(sv1,sv2,sv_names)

    # Plot state-vector comparison
    plot_state_variable_comparison(t,sv1,sv2,sv_names)


if __name__ == "__main__":
    main()
