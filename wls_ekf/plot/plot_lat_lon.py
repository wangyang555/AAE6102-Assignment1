from gnsscommon import *
import numpy as np
import matplotlib.pyplot as plt

def read_filename(filename):
    lat = []
    lon = []
    with open(filename) as fopen:
        lines = fopen.readlines()
        for line in lines:
            line_data = line.split()
            xyz = np.array([float(line_data[1]), float(line_data[2]), float(line_data[3])])
            pos = ecef2pos(xyz)
            lat.append(np.rad2deg(pos[0]))
            lon.append(np.rad2deg(pos[1]))
    return lat, lon


def plot_lat_lon(i,truth_lat, truth_lon, lat, lon, size, title):
    plt.figure(i)
    plt.scatter(lon, lat, label='Data Points', color='blue')
    plt.scatter(truth_lon, truth_lat, label='True Value', color='red', marker='^', s=200)

    plt.xlim(truth_lon - size, truth_lon + size)
    plt.ylim(truth_lat - size, truth_lat + size)
    plt.title(title)
    plt.xlabel('Longitude (degree)')
    plt.ylabel('Latitude (degree)')
    plt.grid(True)

"""
Opensky
22.328444770087565, 114.1713630049711
"""

truth_lat = 22.328444770087565
truth_lon = 114.1713630049711
filename = r"..\result\Opensky\wls.txt"
lat, lon = read_filename(filename)
title = 'Opensky WLS Positioning Result'
plot_lat_lon(1, truth_lat, truth_lon, lat, lon, 0.0001, title)

filename = r"..\result\Opensky\ekf.txt"
lat, lon = read_filename(filename)
title = 'Opensky EKF Positioning Result'
plot_lat_lon(2, truth_lat, truth_lon, lat, lon, 0.0001, title)

"""
Urban
22.3198722, 114.209101777778
"""
truth_lat = 22.3198722
truth_lon = 114.209101777778
filename = r"..\result\Urban\wls.txt"
lat, lon = read_filename(filename)
title = 'Urban WLS Positioning Result'
plot_lat_lon(3, truth_lat, truth_lon, lat, lon, 0.001, title)

filename = r"..\result\Urban\ekf.txt"
lat, lon = read_filename(filename)
title = 'Urban EKF Positioning Result'
plot_lat_lon(4, truth_lat, truth_lon, lat, lon, 0.001, title)

plt.show()
