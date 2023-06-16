#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import pymatreader
import dataPatients

# from os.path import dirname, join as pjoin
# import scipy.io as sio

#  For the induction phase:
#  • TT: observed time-to-target (in seconds) required for reaching
#  the first time the target interval of 45–55 BIS values;
#  • BIS-NADIR: the lowest observed BIS value;
#  • ST10: settling time, defined as the time interval for the BIS to
#  reach and steady within the BIS range between 45 and 55 (that
#  is, the target value of 50 ± 5);
#  • ST20: the same of ST10 but it considers a BIS range of 40 and 60;
#  • US: undershoot, defined as the difference between the lower
#  threshold of 45 and the minimum value of BIS below this threshold.

#  Ffor maintencace phase:
#  • TT: observed time-to-target (in seconds) required for reaching
#  the first time the target interval of 45–55 BIS values;
#  • BIS-NADIR: the lowest observed BIS value;


TT = np.zeros(13)
BIS_NADIR = np.zeros(13)
ST10 = np.zeros(13)
ST20 = np.zeros(13)
US = np.zeros(13)
TTp = np.zeros(13)
BIS_NADIRp = np.zeros(13)
TTn = np.zeros(13)
BIS_NADIRn = np.zeros(13)
IAE = np.zeros(13)

# load simualtion data
mat_file = "data_gpc_dist_profile_2.mat"
sim_data = pymatreader.read_mat(mat_file)

# simulation data following the control scenario from Orignal Paper

Tsim = 20*60
SP = 50  # setpoint value
tstep0 = 0  # setpoint step change time instant
astep = 10  # amplitude of distrubance
tstep1 = 11*60  # possitive disturbance step change time instant
tstep2 = 16*60  # negative disturbance step change time instant

y = sim_data['y']  # BIS signal
u = sim_data['u']  # control signal
t = sim_data['t']  # time vector


for i in range(len(y)):
    istep0 = t[i][tstep0]
    istep1 = t[i][tstep1]
    istep2 = t[i][tstep2]

    BIS_NADIR[i] = min(y[i][0:istep1])
    BIS_NADIRp[i] = min(y[i][istep1:istep2-1])
    BIS_NADIRn[i] = max(y[i][istep2:len(y[i])])

    if BIS_NADIR[i] < 45:
        US[i] = 45-BIS_NADIR[i]
    else:
        US[i] = 0

    flag = 0

    for k in range(istep1-1):
        # if k>70:
        #     stop=1
        if (y[i][k] <= 55) & (flag == 0):
            TT[i] = (t[i][k]-t[i][istep0]) / 60
            flag = 1
        if np.abs(y[i][k]-50) > 5:
            ST10[i] = (t[i][k+1]-t[i][istep0]) / 60
        if np.abs(y[i][k]-50) >= 10:
            ST20[i] = (t[i][k+1]-t[i][istep0]) / 60

    flag = 0
    for k in range(istep1, istep2 - 1, 1):
        if (y[i][k] < 55) & (flag == 0):
            TTp[i] = (t[i][k+1] - t[i][istep1]) / 60
            flag = 1

    flag = 0
    for k in range(istep2, len(y[i]), 1):
        if (y[i][k] > 45) & (flag == 0):
            TTn[i] = (t[i][k+1] - t[i][istep2]) / 60
            flag = 1

    IAE[i] = np.trapz(np.abs(y[i] - SP), t[i])

fig, axs = plt.subplots(2, 1)
fig.suptitle('Simulation results for 13 virtual patients')
for i in range(len(y)):
    axs[0].plot(t[i][:], y[i][:])
    axs[1].plot(t[i][:], u[i][:])
axs[0].set(ylabel='BIS')
axs[1].set(xlabel='time (s)', ylabel='Infusion')
fig.savefig("test.png")
plt.show()

Patients = dataPatients.clinical_data()
print(Patients)

# performance indexes summary
performance_indexes = {"Patient_ID": list(range(1, 14, 1)),
                       # induction phase
                       "TT": TT,
                       "ST20": ST20,
                       "ST10": ST10,
                       "US": US,
                       # maintenance pahse
                       "TTp": TTp,
                       "BIS_NADIRp": BIS_NADIRp,
                       "TTn": TTn,
                       "BIS_NADIRn": BIS_NADIRn,
                       "IAE": IAE}

print(performance_indexes)
