from wls import *
from ekf import *

def output_data_to_file(data, filname):
    """
    output the pos and vel result to file
    """
    out_str = ""
    for i in range(len(data)):
        for j in range(len(data[i])):
            out_str += "%.3f  " % (data[i][j])
        out_str += "\n"
    with open(filname, 'w') as fopen:
        fopen.writelines(out_str)

def wls_method(obs_filename, eph_filename, output_filename):
    """
    wls
    """
    out = WLS_pos_vel_estimation(obs_filename,eph_filename)
    output_data_to_file(out, output_filename)

def ekf_method(obs_filename, eph_filename, output_filename):
    """
    EKF
    """
    out = EKF_pos_vel_estimation(obs_filename, eph_filename)
    output_data_to_file(out, output_filename)

def main():
    """
    outputformat: tow,x,y,z,vx,vy,vz
    """
    # Opensky
    opensky_obs_filename = r".\data\Opensky\obsData.mat"
    opensky_eph_filename = r".\data\Opensky\ephData.mat"
    opensky_output_filename_wls = r".\result\Opensky\wls.txt"
    opensky_output_filename_ekf = r".\result\Opensky\ekf.txt"
    wls_method(opensky_obs_filename, opensky_eph_filename, opensky_output_filename_wls)
    ekf_method(opensky_obs_filename, opensky_eph_filename, opensky_output_filename_ekf)
    """
    # Urban doppler need to be corrected to -doppler
    urban_obs_filename = r".\data\Urban\obsData.mat"
    urban_eph_filename = r".\data\Urban\ephData.mat"
    urban_output_filename_wls = r".\result\Urban\wls.txt"
    urban_output_filename_ekf = r".\result\Urban\ekf.txt"
    wls_method(urban_obs_filename, urban_eph_filename, urban_output_filename_wls)
    ekf_method(urban_obs_filename, urban_eph_filename, urban_output_filename_ekf)
    """

if __name__ == '__main__':
    main()
