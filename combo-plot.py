from RF import RFpower
from DoubleProbe import DoubleProbe
from OceanSpectra import OceanSpectra

import numpy as np
import matplotlib.pyplot as plt
import sys

from glob import glob
from matplotlib.gridspec import GridSpec

'''
This script combines data from NIDAQ, RF, and Spectroscopy (if available)
usage: python comboPlot [shotnumber]

23 December 2023
'''

shot = sys.argv[1]
path = f"data/{shot}/"

hasSpec = True

fig,axs = plt.subplots(5,1, figsize=(12,10))

# double probe
probe = DoubleProbe(path+"NIDAQtext.txt")
probe.plotIV()

# RF power
rf1 = RFpower(path+"RFLog1.txt")
rf1.plotRF()

#rf2 = RFpower(path+"RFLog2.txt")
#rf1.comboPlot(rf2)
T_rf1_fwd = rf1.t_fwd_abs
T_rf1_rev = rf1.t_rev_abs
p1_fwd = rf1.P_fwd
p1_rev = rf1.P_rev

# seconds between Jan 1 1904 and Jan 1 1970, both GMT midnight
t_gap = 2082844800.0


# spectroscopy
if hasSpec:
    path2 = glob(f"data/spectroscopy/*{shot}*txt")[0]
    spec = OceanSpectra(path2)
    spec.plot2d(j=150)
    t_spec = spec.unix_time # ms from 1970
    T_spec = t_spec/1e3  + t_gap

# get common time
t_probe = probe.unix_time # s from 1904


T_probe = t_probe 


if hasSpec:
    t0_global = np.min([T_spec[0], T_probe[0], T_rf1_fwd[0], T_rf1_rev[0]])
    t1_global = np.max([T_spec[-1], T_probe[-1], T_rf1_fwd[-1], T_rf1_rev[-1]]) - t0_global

else:
    t0_global = np.min([T_probe[0], T_rf1_fwd[0], T_rf1_rev[0]])
    t1_global = np.max([T_probe[-1], T_rf1_fwd[-1], T_rf1_rev[-1]]) - t0_global

if hasSpec:
    T_spec -= t0_global
    spec.findPeak(f0=656.279) #H-alpha
    spec.findPeak(f0=486.135) #H-beta
    spec.findPeak(f0=434.0462) #H-gamma

    # plot identified freq peaks over time
    N_lines = len(spec.lines)
    for j in np.arange(N_lines):
        axs[0].plot(T_spec, spec.lines[j], label=f"{spec.freqs[j]} nm")
    axs[0].set_ylabel('counts')

axs[0].set_title(shot)

T_probe -= t0_global
T_rf1_fwd -= t0_global
T_rf1_rev -= t0_global


probe.plotPressure(axs[2])

# rf data
axs[1].plot(T_rf1_fwd,p1_fwd,'o-',label="P forward")
axs[1].plot(T_rf1_rev,p1_rev,'o-',label="P reverse")

axs[3].plot(T_probe, probe.I, label='mA')
axs[4].plot(T_probe, probe.V, label='mA')

axs[-1].set_xlabel('time (s)')
for a in axs:
    a.set_xlim(0,t1_global) 
    a.legend()
    a.grid()

fig.tight_layout()
plt.show()
