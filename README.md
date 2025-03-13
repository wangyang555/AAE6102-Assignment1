# Satellite Communication and Navigation Assignment (AAE6102)

**Name:** WANG Yang  
**Email:** 24041885r@connect.polyu.hk  

This repository contains the implementation of Assignment 1 for AAE6102 Satellite Communication and Navigation. The assignment consists of 5 tasks. The first three tasks were completed using the open-source MATLAB project [FGI-GSRx](https://github.com/nlsfi/FGI-GSRx). For Tasks 4-5, I implemented the solutions in Python, utilizing some functions from the [rtklibexplorer/rtklib-py](https://github.com/rtklibexplorer/rtklib-py), such as ionospheric delay, tropospheric delay, and satellite position calculations. Based on these functions, I developed the **Weighted Least Squares (WLS)** and **Extended Kalman Filter (EKF)** algorithms to estimate the receiver's position and velocity using pseudorange and Doppler measurements. Feedback and discussions are welcome!

---

## Task 1 — Acquisition  
Process IF data using a GNSS SDR and generate initial acquisition outputs.  

The signal information of Opensky and Urban is as follows：
| Parameter              | Opensky           | Urban             |  
|------------------------|-------------------|-------------------|  
| **Carrier Frequency**  | 1575.42 MHz       | 1575.42 MHz       |  
| **Intermediate Freq.** | 4.58 MHz          | 0                 |  
| **Sampling Freq.**     | 58 MHz            | 26 MHz            |  
| **Data Format**        | 8-bit I/Q samples | 8-bit I/Q samples |  
| **Ground Truth**       | (22.32844477, 114.17136300) | (22.3198722, 114.20910178) |  
| **Data Length**        | 90 seconds        | 90 seconds        |  

### Results  
- **Opensky Signal acquisition results**
<div align="center">
  <img src="figure/Opensky/1_1_freq_codephase_signal_search.png" alt="">
</div>
<div align="center">
  <img src="figure/Opensky/1_2_codephase_signal_search.png" alt="">
</div>
<div align="center">
  <img src="figure/Opensky/1_3_acquisition_result.png" alt="">
</div>

- **Urban Signal acquisition results**
<div align="center">
  <img src="figure/Urban/1_1_freq_codephase_signal_search.png" alt="">
</div>
<div align="center">
  <img src="figure/Urban/1_2_codephase_signal_search.png" alt="">
</div>
<div align="center">
  <img src="figure/Urban/1_3_acquisition_result.png" alt="">
</div>

---

## Task 2 — Tracking  
Adapt the tracking loop (DLL) to produce correlation plots and analyze tracking performance under urban interference.  

### Results
- **Opensky Signal tracking results**
<div align="center">
  <img src="figure/Opensky/2_1_tracking.png" alt="">
</div>
<div align="center">
  <img src="figure/Opensky/2_2_correlation.png" alt="">
</div>

- **Urban Signal tracking results**
<div align="center">
  <img src="figure/Urban/2_1_tracking.png" alt="">
</div>
<div align="center">
  <img src="figure/Urban/2_2_correlation.png" alt="">
</div>

### Analysis:  
- **Opensky:** Correlation peaks exhibit a symmetrical triangular shape.  
- **Urban:** Multipath effects and NLOS signals cause correlation peak distortions (e.g., multiple peaks, shifted/broadened main peaks) and affect the measurement accuracy of the code phase, reducing code-phase measurement accuracy.  

---

## Task 3 — Navigation Data Decoding  
Decode navigation messages and extract ephemeris data.  

### Results
- **Opensky Navigation messages decoded**
Navigation messages decoded for **9 satellites**
<div align="center">
  <img src="figure/Opensky/3_1_eph_message.png" alt="">
</div>
<div align="center">
  <img src="figure/Opensky/3_2_cn0.png" alt="">
</div>

- **Urban Navigation messages decoded**
Navigation messages decoded for **4 satellites** (lower C/N₀ values due to urban interference)
<div align="center">
  <img src="figure/Urban/3_1_eph.png" alt="">
</div>
<div align="center">
  <img src="figure/Urban/3_2_cn0.png" alt="">
</div>

---

## Task 4 — Position & Velocity Estimation (WLS)  
Implement WLS using pseudorange measurements.  

### Model:
Below I give the function model and Stochastic model of pseudorange and doppler observations. For more information, please refer to the file under document.
<div align="center">
  <img src="document/img_for_github/Pseudorange.png" alt="">
</div>
<div align="center">
  <img src="document/img_for_github/doppler.png" alt="">
</div>

### Results:  
1. **Position Comparison:** Longitude/latitude vs. ground truth.

**Opensky**
<div align="center">
  <img src="figure/Opensky/4_1_opensky_wls_lat_lon.png" alt="">
</div>

**Urban**
<div align="center">
  <img src="figure/Urban/4_1_urban_wls_lat_lon.png" alt="">
</div>

2. **NEU Time Series:** Relative to epoch-averaged coordinates.

**Opensky**
<div align="center">
  <img src="figure/Opensky/4_2_opensky_wls_neu.png" alt="">
</div>

**Urban**
<div align="center">
  <img src="figure/Urban/4_2_urban_wls_neu.png" alt="">
</div>

3. **XYZ Velocity Time Series.**  

**Opensky**
<div align="center">
  <img src="figure/Opensky/4_3_opensky_wls_v.png" alt="">
</div>

**Urban**
<div align="center">
  <img src="figure/Urban/4_3_urban_wls_v.png" alt="">
</div>

### Analysis:  
It can be seen that in the Urban environment, the accuracy is significantly lower than that in the Opensky environment. On the one hand, there are fewer satellites tracked in the Urban environment. On the other hand, the pseudorange does not take into account the impact of multipath. Complex reflective surfaces in the Urban environment, such as glass curtain walls, cause serious multipath interference, which has a more obvious negative impact on positioning.

---

## Task 5 — Kalman Filter-Based Positioning (EKF)  
Develop an EKF using pseudorange and Doppler measurements.  

### Results:  
1. **Position Comparison:** Longitude/latitude vs. ground truth.
**Opensky**
<div align="center">
  <img src="figure/Opensky/5_1_opensky_ekf_lat_lon.png" alt="">
</div>

**Urban**
<div align="center">
  <img src="figure/Urban/5_1_urban_ekf_lat_lon.png" alt="">
</div>

2. **NEU Time Series:** Relative to epoch-averaged coordinates.
**Opensky**
<div align="center">
  <img src="figure/Opensky/5_2_opensky_ekf_neu.png" alt="">
</div>

**Urban**
<div align="center">
  <img src="figure/Urban/5_2_urban_ekf_neu.png" alt="">
</div>

3. **XYZ Velocity Time Series.**  
**Opensky**
<div align="center">
  <img src="figure/Opensky/5_3_opensky_ekf_v.png" alt="">
</div>

**Urban**
<div align="center">
  <img src="figure/Urban/5_3_urban_ekf_v.png" alt="">
</div>

---

## References  
1. [FGI-GSRx](https://github.com/nlsfi/FGI-GSRx](https://github.com/nlsfi/FGI-GSRx)  
2. [rtklibexplorer/rtklib-py](https://github.com/rtklibexplorer/rtklib-py](https://github.com/rtklibexplorer/rtklib-py)  

---

## Note  
The EKF was partially assisted by [**DeepSeek**](https://chat.deepseek.com/), an AI tool, for code optimization and debugging.  


