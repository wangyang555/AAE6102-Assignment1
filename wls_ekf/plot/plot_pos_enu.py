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
            xyz.append([float(line_data[1]), float(line_data[2]), float(line_data[3])])
    return xyz

def calculate_average(xyz):
    xyz_array = np.array(xyz)
    x_mean = np.mean(xyz_array[:, 0])
    y_mean = np.mean(xyz_array[:, 1])
    z_mean = np.mean(xyz_array[:, 2])
    return x_mean, y_mean, z_mean

def calculate_neu(x_mean,y_mean,z_mean,xyz):
    neu = []
    for i in range(len(xyz)):
        n, e, u = xyz2neu(x_mean, y_mean, z_mean, xyz[i][0], xyz[i][1], xyz[i][2])
        neu.append([n,e,u])
    return neu


def plot_lat_lon(i,neu, size, title):
    neu = np.array(neu)
    n = neu[:, 0]
    e = neu[:, 1]
    u = neu[:, 2]
    indices = range(len(n))
    plt.figure(i,figsize=(10, 6))
    plt.ylim(-size, size)
    plt.plot(indices, n, marker='o', linestyle='-', color='r', label='N Coordinates (m)')
    plt.plot(indices, e, marker='s', linestyle='-', color='g', label='E Coordinates (m)')
    plt.plot(indices, u, marker='^', linestyle='-', color='b', label='U Coordinates (m)')

    plt.legend()
    plt.xlabel('Index')
    plt.ylabel('Coordinate Value')
    plt.title(title)
    plt.grid(True)



# Opensky
filename = r"..\result\Opensky\wls.txt"
xyz = read_filename(filename)
x_mean, y_mean, z_mean = calculate_average(xyz)
neu = calculate_neu(x_mean,y_mean,z_mean,xyz)
title = 'Opensky WLS NEU Result'
plot_lat_lon(1, neu, 20,title)

filename = r"..\result\Opensky\ekf.txt"
xyz = read_filename(filename)
x_mean, y_mean, z_mean = -2416982.251,5385230.181,2408088.109 # The first epoch influence the average value
neu = calculate_neu(x_mean,y_mean,z_mean,xyz)
title = 'Opensky EKF NEU Result'
plot_lat_lon(2, neu, 20, title)


# Urban
filename = r"..\result\Urban\wls.txt"
xyz = read_filename(filename)
x_mean, y_mean, z_mean = calculate_average(xyz)
neu = calculate_neu(x_mean,y_mean,z_mean,xyz)
title = 'Urban WLS NEU Result'
plot_lat_lon(3, neu, 300, title)

filename = r"..\result\Urban\ekf.txt"
xyz = read_filename(filename)
x_mean, y_mean, z_mean = -2420691.890,5384037.182,2407276.987
neu = calculate_neu(x_mean,y_mean,z_mean,xyz)
title = 'Urban EKF NEU Result'
plot_lat_lon(4, neu, 300, title)

plt.show()