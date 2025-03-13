from gnsscommon import *
import numpy as np
import matplotlib.pyplot as plt
from pycoord import *

def read_filename(filename):
    xyz = []
    with open(filename) as fopen:
        lines = fopen.readlines()
        for line in lines:
            line_data = line.split()
            xyz.append([float(line_data[4]), float(line_data[5]), float(line_data[6])])
    return xyz


def plot_lat_lon(i,xyz, size, title):
    xyz = np.array(xyz)
    x = xyz[:, 0]
    y = xyz[:, 1]
    z = xyz[:, 2]
    indices = range(len(x))
    plt.figure(i,figsize=(10, 6))
    plt.ylim(-size, size)
    plt.plot(indices, x, marker='o', linestyle='-', color='r', label='X Velocity (m/s)')
    plt.plot(indices, y, marker='s', linestyle='-', color='g', label='Y Velocity (m/s)')
    plt.plot(indices, z, marker='^', linestyle='-', color='b', label='Z Velocity (m/s)')

    plt.legend()
    plt.xlabel('Index')
    plt.ylabel('Coordinate Value')
    plt.title(title)
    plt.grid(True)



# Opensky
filename = r"..\result\Opensky\wls.txt"
xyz = read_filename(filename)
title = 'Opensky WLS XYZ Velocity'
plot_lat_lon(1, xyz, 2,title)

filename = r"..\result\Opensky\ekf.txt"
xyz = read_filename(filename)
title = 'Opensky EKF XYZ Velocity'
plot_lat_lon(2, xyz, 2, title)


# Urban
filename = r"..\result\Urban\wls.txt"
xyz = read_filename(filename)
title = 'Urban WLS XYZ Velocity'
plot_lat_lon(3, xyz, 3000, title)

filename = r"..\result\Urban\ekf.txt"
xyz = read_filename(filename)
title = 'Urban EKF XYZ Velocity'
plot_lat_lon(4, xyz, 3000, title)

plt.show()