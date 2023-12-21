import numpy as np
import matplotlib.pyplot as plt
import sys

from matplotlib.gridspec import GridSpec
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit

'''
Muse Data Analysis Script

T Qian, 19 Dec 2023
'''


class DoubleProbe():

    def __init__(self,fin):

        self.fname = fin
        self.loadData(fin)


    def loadData(self,fin):

        data = np.genfromtxt(fin,delimiter=',')

        time_long, AI0, AI1, AI2, AI3, AO0, AO1 = data.T  # s
        self.time_unix = time_long
        self.AI0 = AI0  # pressure
        self.AI1 = AI1  # bias request
        self.AI2 = AI2  # bias
        self.AI3 = AI3  # shunt
        self.AO0 = AO0  # bias, same as AI2
        self.AO1 = AO1  # unused

        self.time = time_long - time_long[0] # seconds

    def filterIV(self, t0=2, # crop start time
                            t1=7, # crop finish time
                            sg_window = 50, # savgol window
                            sg_order = 3, # savgol polynomial order
                      ):

        # should code the conversion in ONE place to avoid bugs
        V_probe = self.AI2*10  # V
        I_probe = self.AI3/5*97.8757 # mA
        time = self.time # s

        def argNear(arr,val):
            return np.argmin( np.abs(arr - val) )

        t0_idx = argNear(time,t0)
        t1_idx = argNear(time,t1)

        # cut
        V_cut = V_probe[t0_idx:t1_idx]
        I_cut = I_probe[t0_idx:t1_idx]
        t_cut = time[t0_idx:t1_idx]

        # filter
        V_filter = savgol_filter(V_cut, sg_window, sg_order)
        I_filter = savgol_filter(I_cut, sg_window, sg_order)

        # fit
        def IV_tanh(Vbias,Te,Isat,I_offset):
            return Isat * np.tanh( Vbias / 2 / Te ) + I_offset
        param, cov = curve_fit(IV_tanh, V_filter, I_filter)
        param2, cov2 = curve_fit(IV_tanh, V_cut, I_cut)

        I_fit = IV_tanh(V_filter,*param)
        I_fit2 = IV_tanh(V_filter,*param2)


        # plot
        axs = self.ax_IV
        axs[0].plot(t_cut, V_cut, 'C3.', ms=1)
        axs[1].plot(t_cut, I_cut, 'C3.', ms=1)
        axs[2].plot(V_cut, I_cut, 'C3.', ms=1)

        axs[0].plot(t_cut, V_filter, 'C2--', lw=0.5)
        axs[1].plot(t_cut, I_filter, 'C2--', lw=0.5)
        axs[2].plot(V_cut, I_filter, 'C2--', lw=2, label=f"(win,order) = {sg_window}, {sg_order}")

        Te, Isat, I_offset = param
        axs[2].plot(V_cut, I_fit, 'C1:', ms=12, label=f"(Te, Isat, Ioffset) = {Te:.1f} eV, {Isat:.3f} mA, {I_offset:.2f} mA")

        Te, Isat, I_offset = param2

        err2 = np.sqrt( np.diag(cov2) )
        dT, dIsat, dIoff = err2
        axs[2].plot(V_cut, I_fit2, 'C1--', ms=12, label="fit")
        axs[2].plot([],[],' ', label=rf"Te = {Te:.1f}$\pm${dT:.1f} eV")
        axs[2].plot([],[],' ', label=rf"Isat = {Isat:.2f}$\pm${dIsat:.2f} mA")
        axs[2].plot([],[],' ', label=rf"I_offset = {I_offset:.2f}$\pm${dIoff:.2f} mA")

        axs[2].legend()

        axs[2].set_ylim( 1.1*np.min(I_cut), 1.1*np.max(I_cut) )


    def plotRaw(self, save=False):

        time_long = self.time_unix
        pressure = self.AI0
        V_bias_request = self.AI1 
        V_bias = self.AI2
        V_shunt = self.AI3


        # raw data 
        fig,axs = plt.subplots(5,1, figsize=(10,6))
        axs[0].plot(time_long, label='time')
        axs[1].plot(pressure, label='pressure')
        axs[2].plot(V_bias_request, label='bias request') #5x
        axs[3].plot(V_bias, label='bias') # 10x
        axs[4].plot(V_shunt, label='shunt') # divide by 5, multiply 97 gain, divide by 1kOhm
        [a.legend(loc=0) for a in axs]
        fig.suptitle(self.fname)
        plt.tight_layout()

        if save:
            fig.savefig(save)


    def plotIV(self, save=False):

        fig = plt.figure(layout="constrained", figsize=(10,6))
        gs = GridSpec(2, 2, figure=fig)
        ax1 = fig.add_subplot(gs[0, 1])
        ax2 = fig.add_subplot(gs[1, 1])
        ax3 = fig.add_subplot(gs[:, 0])
        
        V_probe = self.AI2*10  # V
        I_probe = self.AI3/5*97.8757 # mA
        time = self.time # s
        ax1.plot(time, V_probe, 'C0.', label='V')
        ax2.plot(time, I_probe, 'C1.', label='mA')
        ax3.plot(V_probe, I_probe, 'C4.')

#        if np.min(I_probe) < -5:
#            m = np.max(I_probe)
#            ax3.set_ylim(-1, 1.1*m)
        ax2.set_xlabel('time [s]')
        ax3.set_ylabel('mA')
        ax3.set_xlabel('V')
        [a.legend(loc=0) for a in [ax1,ax2] ]
        [a.grid() for a in [ax1,ax2,ax3] ]
                              
        fig.suptitle(self.fname)
        fig.tight_layout()

        if save:
            fig.savefig(save)

        self.ax_IV = [ax1,ax2,ax3]

    def plotPressure(self, save=False):

        V = self.AI0
        time = self.time

        P_raw = 10**((V - 5.5)/0.5)
        P_H2 = P_raw / 0.42

        fig,axs = plt.subplots(1,1)
        axs.plot(time, P_raw, '.', label="pressure raw")
        axs.plot(time, P_H2, '.', label="pressure H2")
  
        axs.set_ylabel('Torr')
        axs.set_xlabel('s')
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0) )
        # how to make y-axis in sci notation?
        # plt.gca().ticklabel_format(axis='y', style='sci',


        plt.legend()
        plt.grid()

        fig.suptitle(self.fname)
        fig.tight_layout()

        if save:
            fig.savefig(save)


if __name__ == '__main__':
    path = sys.argv[1]
    data = DoubleProbe(path+"NIDAQtext.txt")
    
    data.plotRaw()
    data.plotIV()
    data.filterIV()
    data.plotPressure()
    
    plt.show()
    

