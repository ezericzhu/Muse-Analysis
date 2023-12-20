import numpy as np
import matplotlib.pyplot as plt
import sys

from matplotlib.gridspec import GridSpec

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

        ax2.set_xlabel('time [s]')
        ax3.set_ylabel('mA')
        ax3.set_xlabel('V')
        [a.legend(loc=0) for a in [ax1,ax2] ]
        [a.grid() for a in [ax1,ax2,ax3] ]
                              
        fig.suptitle(self.fname)
        fig.tight_layout()

        if save:
            fig.savefig(save)

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
    data.plotPressure()
    
    plt.show()
    

