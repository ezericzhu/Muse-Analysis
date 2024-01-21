from .RF import RFpower
from .DoubleProbe import DoubleProbe
from .OceanSpectra import OceanSpectra

import numpy as np
import matplotlib.pyplot as plt
# import sys
import os

from glob import glob
from matplotlib.gridspec import GridSpec
import argparse

def main():
    '''
    This script combines data from NIDAQ, RF, and Spectroscopy (if available)

    usage: python MuseAnalysis/comboPlot.py [shotnumber]
    Note paths, default assumes script is called from folder with data/231223000

    23 December 2023
    '''

    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument("-t","--target")
    parser.add_argument("-s","--show",action="store_true")

    args = parser.parse_args()

    if hasattr(args, "help"):
        if args.help:
            parser.print_help()
    if hasattr(args, "target"):
        shot = args.target
    else:
        raise ValueError("No target specified")
    if hasattr(args, "show"):
        iflag_show = args.show
    else:
        iflag_show = False

    # shot = sys.argv[1]
    
    data_path = os.getenv("MUSE_DATA_PATH")
    if data_path is None:
        data_path = "./"
    path = f"{data_path}{shot}/"
    

    fig,axs = plt.subplots(5,1, figsize=(12,9))
    axs[0].set_title(shot)

    # seconds between Jan 1 1904 and Jan 1 1970, both GMT midnight
    t_gap = 2082844800.0

    ### Load Data

    # double probe
    try:
        probe = DoubleProbe(path+"NIDAQtext.txt")

        if int(shot) < 231222000: # threshold for bias box change
            probe.V_factor = 10
            probe.I_factor = 97.8757/5

        probe.plotIV(save=path+"plotIV.png")
        probe.plotRaw(save=path+"plotRaw.png")

    except:
        print("no probe data")

    # spectroscopy
    try:
        path2 = glob(f"data/spectroscopy/*{shot}*txt")[0]
        spec = OceanSpectra(path2)
        t_spec = spec.unix_time # ms from 1970
        T_spec = t_spec/1e3  + t_gap # sec from 1904
        hasSpec = True
    except:
        print("spectroscopy data not found")
        hasSpec = False

    # RF power
    try:
        rf1 = RFpower(path+"RFLog1.txt")
        rf1.plotRF()
        hasRF1 = True
    except:
        hasRF1 = False

    try:
        rf2 = RFpower(path+"RFLog2.txt")
        rf2.plotRF()
        hasRF2 = True
    except:
        hasRF2 = False

    if hasRF1 and hasRF2:
        rf1.comboPlot(rf2)



    ### Get Common Time
    T_probe = probe.unix_time # s from 1904


    # this will break if there is no RF1 data
    T_rf1_fwd = rf1.t_fwd_abs
    T_rf1_rev = rf1.t_rev_abs
    p1_fwd = rf1.P_fwd
    p1_rev = rf1.P_rev
    t0_global = T_rf1_fwd[0]

    T_probe -= t0_global
    T_rf1_fwd -= t0_global
    T_rf1_rev -= t0_global
    if hasSpec:
        T_spec -= t0_global

    if hasSpec:
        t_start = np.min([T_spec[0], T_probe[0], T_rf1_fwd[0], T_rf1_rev[0]])
        t_end = np.max([T_spec[-1], T_probe[-1], T_rf1_fwd[-1], T_rf1_rev[-1]]) #- t0_global
    else:
        t_start = np.min([T_probe[0], T_rf1_fwd[0], T_rf1_rev[0]])
        t_end = np.max([T_probe[-1], T_rf1_fwd[-1], T_rf1_rev[-1]]) #- t0_global


    ### Plot

    if hasSpec:
        spec.findPeak(f0=656.279) #H-alpha
        spec.findPeak(f0=486.135) #H-beta
        spec.findPeak(f0=434.0462) #H-gamma

        # plot identified freq peaks over time
        N_lines = len(spec.lines)
        for j in np.arange(N_lines):
            axs[0].plot(T_spec, spec.lines[j], label=f"{spec.freqs[j]} nm")
        axs[0].set_ylabel('counts')

        # make a second panel with 2D spectragram
        spec.plot2d(j=150)


    axs[1].plot(T_rf1_fwd, p1_fwd,'o-',label="P forward")
    axs[1].plot(T_rf1_rev, p1_rev,'o-',label="P reverse")

    probe.plotPressure(axs[2], t_global=T_probe)
    # probe.plotPressure()
    axs[2].grid()

    axs[3].plot(T_probe, probe.V, label='probe V')
    axs[4].plot(T_probe, probe.I, label='shunt I')

    axs[1].set_ylabel('RF Power (W)')
    axs[3].set_ylabel('V')
    axs[4].set_ylabel('mA')

    for a in axs:
        a.set_xlim(t_start, t_end) 
        a.legend()
        a.grid()
    axs[-1].set_xlabel('time (s)')

    fig.tight_layout()
    fig.savefig(path+"plotTime.png")
    if iflag_show:
        plt.show()


if __name__ == "__main__":
    main()
