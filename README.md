# Satellite Communication and Navigation Assignment (AAE6102)

**Name:** WANG Yang  
**Email:** 24041885r@connect.polyu.hk  

This repository contains the implementation of **Assignment 1** for AAE6102 Satellite Communication and Navigation. The assignment includes 5 tasks, with the first three tasks processed using the open-source MATLAB project [FGI-GSRx](https://github.com/nlsfi/FGI-GSRx), and Tasks 4-5 implemented in Python. For Tasks 4-5, functions from the [rtklib-py](https://github.com/rtklibexplorer/rtklib-py) repository (e.g., ionospheric/tropospheric delay calculations) were utilized to develop **Weighted Least Squares (WLS)** and **Extended Kalman Filter (EKF)** algorithms for receiver position/velocity estimation. Feedback and discussions are welcome!

---

## Task 1 — Acquisition  
Process IF data using a GNSS SDR and generate initial acquisition outputs.  

The signal information of Opensky and Urban is as follows
| Parameter              | Opensky           | Urban             |  
|------------------------|-------------------|-------------------|  
| **Carrier Frequency**  | 1575.42 MHz       | 1575.42 MHz       |  
| **Intermediate Freq.** | 4.58 MHz          | 0                 |  
| **Sampling Freq.**     | 58 MHz            | 26 MHz            |  
| **Data Format**        | 8-bit I/Q samples | 8-bit I/Q samples |  
| **Ground Truth**       | (22.32844477, 114.17136300) | (22.3198722, 114.20910178) |  
| **Data Length**        | 90 seconds        | 90 seconds        |  

**Results:**  
- **Opensky:** Signal acquisition results (see Figures 1-3).  
- **Urban:** Signal acquisition results (see Figures 1-3).  

---

## Task 2 — Tracking  
Adapt the tracking loop (DLL) to produce correlation plots and analyze tracking performance under urban interference.  

### Key Observations:  
- **Opensky:** Correlation peaks exhibit a symmetrical triangular shape.  
- **Urban:** Multipath effects and NLOS signals cause correlation peak distortions (e.g., multiple peaks, shifted/broadened main peaks), reducing code-phase measurement accuracy.  

---

## Task 3 — Navigation Data Decoding  
Decode navigation messages and extract ephemeris data.  

### Results:  
- **Opensky:** Navigation messages decoded for **9 satellites**.  
- **Urban:** Navigation messages decoded for **4 satellites** (lower C/N₀ values due to urban interference).  

---

## Task 4 — Position & Velocity Estimation (WLS)  
Implement WLS using pseudorange measurements.  

### Results:  
#### Opensky Environment  
1. **Position Comparison:** Longitude/latitude vs. ground truth.  
2. **NEU Time Series:** Relative to epoch-averaged coordinates.  
3. **XYZ Velocity Time Series.**  

#### Urban Environment  
- **Lower Accuracy:** Fewer tracked satellites and unmitigated multipath effects (e.g., reflections from glass surfaces).  

---

## Task 5 — Kalman Filter-Based Positioning (EKF)  
Develop an EKF using pseudorange and Doppler measurements.  

### Results:  
- **Opensky & Urban:** Position/velocity estimates provided (similar structure to Task 4).  

---

## Notes  
- Figures referenced in the original document (e.g., "Figure.1") are placeholders. Actual figures will be added to the repository.  
- Code implementations and configuration files will be uploaded separately.  

For detailed equations and methodology, refer to the documentation in the repository.  
