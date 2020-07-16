import sys
import numpy as np

filename = sys.argv[1]
data = np.genfromtxt(filename)

mean = np.mean(data)
std = np.std(data)

print("%g +/- %g" % (mean,std))
