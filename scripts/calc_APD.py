from sys import argv

def forwarddiff(y, h):
    n = len(y)
    res = []
    i = 1

    for i in range(1,n):
        res.append((y[i] - y[i-1]) / h);

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

    d = data[start:(start+end)]

    max_v = max(d)

    max_index = d.index(max_v) + start

    return max_index

def index_activation(data, start=0):
    d = data[start:]

    for i, v in enumerate(d):
        if d[i+start] < 0.0 and d[i+start+1] > 0.0:
            return i+start

ap_file_name = argv[1]
num_aps = int(argv[2])
ms_each_step = float(argv[3])

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

for j in range(num_aps):

    slope_s = slope_start(ap_data, slope_e, 0.1)
    max_ind = max_index(ap_data, slope_e, 1050)
    slope_e = slope_end(ap_data, max_ind+3, 0.005)

    apds.append((slope_e-slope_s)*ms_each_step)

print("APD: ", sum(apds)/len(apds))

