# ==================================================================================================
# Author: Lucas Berg @bergolho
# Date: 10/11/2023
# Script to combine the ALG and FIBERS file into a single one.
# --------------------------------------------------------------------------------------------------
# Dependencies: Pandas and NumPy
# ==================================================================================================

import sys
import numpy as np
import pandas as pd

def main():
    if len(sys.argv) != 3:
        print("-------------------------------------------------------------------------")
        print("Usage:> python %s <input_alg_file> <input_fibers_file>" % sys.argv[0])
        print("-------------------------------------------------------------------------")
        print("<input_alg_file> = Input ALG file with center, discretization and extra data")
        print("<input_fibers> = Input FIBERS file with (f,s,n)")
        print("-------------------------------------------------------------------------")
        return 1

    input_alg_filename = sys.argv[1]
    input_fibers_filename = sys.argv[2]

    alg_data = pd.read_csv(input_alg_filename, sep=',')
    #print(alg_data)

    fibers_data = pd.read_csv(input_fibers_filename, sep=' ')
    #print(fibers_data)

    result_data = pd.concat([alg_data, fibers_data], axis=1)
    #print(result_data)

    result_data.to_csv("output.alg", header=False, index=False)

if __name__ == "__main__":
	main()
