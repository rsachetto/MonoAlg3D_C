import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
	
    if len(sys.argv) != 2:
        print("------------------------------------------------------------------------------------------------------")
        print("Usage:> python %s <result_filename>" % sys.argv[0])
        print("------------------------------------------------------------------------------------------------------")
        print("Example:> python %s ../outputs/elnaz-purkinje-coupled-errors/results.dat" % (sys,argv[0]))
        print("------------------------------------------------------------------------------------------------------")
        return 1

    input_file = sys.argv[1]

    result_data = np.genfromtxt(input_file)

    data_scenario_1 = result_data[0:20]
    data_scenario_2 = result_data[20:40]

    for i in range(len(data_scenario_1)):
        #print("Individual %2d -- %10g -- %10g" % (i,data_scenario_1[i],data_scenario_2[i]))
        print("%.10lf %.10lf" % (data_scenario_1[i],data_scenario_2[i]))

if __name__ == "__main__":
	main()
