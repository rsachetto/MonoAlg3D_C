import glob
from sys import argv
import numpy as np

dir1 = argv[1]
dir2 = argv[2]

files1  = glob.glob(dir1+"/V_t*")
files2  = glob.glob(dir2+"/V_t*")
data_final = []

for f1, f2 in zip(files1, files2):
    data1 = open(f1, "r")
    data2 = open(f2, "r")
    
    data_list1 = []
    data_list2 = []
    
    for line in data1:
        line = line.strip().split(",")
        data_list1.append(float(line[4]))

    for line in data2:
        line = line.strip().split(",")
        data_list2.append(float(line[4]))

    x =np.array(data_list1)
    y =np.array(data_list2)

    error = (np.sqrt(sum( np.power(x-y,2) )) / np.sqrt( sum( np.power(y,2))))*100.0
    data_final.append(error)


z = np.array(data_final)
print z.mean(axis=None)
