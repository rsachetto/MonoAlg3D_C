import sys
import glob, os
import numpy as np
import matplotlib.pyplot as plt

def parse_file (filename):
    sv = np.genfromtxt(filename)
    n = len(sv)

    print("-------------------------------------------------------------------")
    print("\nFilename :> %s\n" % filename)
    print("-------------------------------------------------------------------")
    print("CPU")
    for i in range(n):
        value = sv[i]
        print("sv[%d] = %g;\n" % (i,sv[i]))
    print("********************************************************************")
    print("GPU")
    for i in range(n):
        value = sv[i]
        print("*((real *) ((char *) sv + pitch * %d) + threadID) = %g;\n" % (i,sv[i]))

def main():

    if len(sys.argv) != 2:
    	print("-------------------------------------------------------------------------")
    	print("Usage:> python %s <input_folder>" % sys.argv[0])
    	print("-------------------------------------------------------------------------")
    	print("<input_folder> = Input folder")
    	print("-------------------------------------------------------------------------")
    	return 1

    input_folder = sys.argv[1]
    os.chdir(input_folder)
    for filename in glob.glob("*.txt"):
        parse_file(filename)

if __name__ == "__main__":
	main()
