import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.gridspec import GridSpec

'''
Muse Data Analysis Script

T Qian, 19 Dec 2023, midnight edition
'''

class RFpower:

    def __init__(self,fin):

        self.fname = fin
        self.loadData(fin)


    def loadData(self,fin):


        with open(fin) as f:
            datain = f.readlines()


        arr = []
        for line in datain:

            data = getPair(line)

            try:
                if len(data) == 2:
                    arr.append(data)
            except:
                continue

        time,power = np.transpose(arr)

        N = len(arr)//2
        t_fwd, t_rev = time.reshape(N,2).T
        P_fwd, P_rev = power.reshape(N,2).T

        # get time offset from first data point
        t0 = float(datain[0][:14])

        # save
        self.T_fwd = t_fwd - t0
        self.P_fwd = P_fwd
        self.T_rev = t_rev - t0
        self.P_rev = P_rev

        self.data = np.array(arr)
        self.t0 = t0

    def plotRF(self):

        fig = plt.figure(layout="constrained")
        gs = GridSpec(3, 1, figure=fig)
        ax0 = fig.add_subplot(gs[:-1])
        ax1 = fig.add_subplot(gs[-1])

        s = self
        ax0.plot(s.T_fwd, s.P_fwd, 'o-', label="Forward Power")
        ax0.plot(s.T_rev, s.P_rev, 'o-', label="Reflected Power")
        ax1.plot(s.T_rev, s.P_rev, 'C1o-', label="Reflected Power")
        
        ax0.set_ylabel("Power (W)")
        ax1.set_ylabel("Reflected (W)")
        ax1.set_xlabel("Time (s)")
        
        ax0.grid()
        ax1.grid()
        ax0.legend()
        
        fig.suptitle(self.fname)
        fig.tight_layout()

        return fig

    def comboPlot(self, rf2):

        rf1 = self
    
        # load data
        t1_fwd = rf1.T_fwd
        t1_rev = rf1.T_rev
        t2_fwd = rf2.T_fwd
        t2_rev = rf2.T_rev
    
        p1_fwd = rf1.P_fwd
        p1_rev = rf1.P_rev
        p2_fwd = rf2.P_fwd
        p2_rev = rf2.P_rev
    
        # plot
        fig = plt.figure(layout="constrained", figsize=(10,8))
        gs = GridSpec(3, 1, figure=fig)
        ax0 = fig.add_subplot(gs[:-1])
        ax1 = fig.add_subplot(gs[-1])
    
        # need to interpolate
        t_axis = np.linspace( rf1.data[0,0], rf1.data[-1,0], 100 ) - rf1.t0
        p1_fwd_interp = np.interp(t_axis, t1_fwd, p1_fwd)
        p2_fwd_interp = np.interp(t_axis, t2_fwd, p2_fwd)
        P_fwd_total = p1_fwd_interp + p2_fwd_interp
        p1_rev_interp = np.interp(t_axis, t1_rev, p1_rev)
        p2_rev_interp = np.interp(t_axis, t2_rev, p2_rev)
        P_rev_total = p1_rev_interp + p2_rev_interp
    
        ax0.plot(t_axis, P_fwd_total, 'C2', lw=3, label="Total Foward")
        ax0.plot(t_axis, P_rev_total, 'C3', lw=3, label="Total Reflected")
        ax0.plot(t1_fwd, p1_fwd, 'C0o--', mfc='none', label="Forward 1")
        ax0.plot(t2_fwd, p2_fwd, 'C0x--', label="Forward 2")
        ax0.plot(t1_rev, p1_rev, 'C1o--', mfc='none', label="Reflected 1")
        ax0.plot(t2_rev, p2_rev, 'C1x--', label="Reflected 2")
    
        ax1.plot(t1_rev, p1_rev, 'C1o--', mfc='none', label="Reflected Power 1")
        ax1.plot(t2_rev, p2_rev, 'C1x--', label="Reflected Power 2")
        
        ax0.set_ylabel("Power (W)")
        ax1.set_ylabel("Reflected (W)")
        ax1.set_xlabel("Time (s)")
        
        ax0.grid()
        ax1.grid()
        ax0.legend()
        
        fig.suptitle(self.fname)
        fig.tight_layout()
    
        return fig
        
##
# helper function
def getPair(line):

    try:
        # if csv
        data = np.array(line.strip().split(','), float)
        return data

    except:


        try:
            # assume time has 14 char
            t = line.strip()[:14]
            p = line.strip()[14:]
            data = np.array([t,p], float)
            return data

        except:

            return False



##

if __name__ == "__main__":
    path = sys.argv[1]
    
    try:
        d1 = RFpower(path+"RFLog1.txt")
        d1.plotRF()
    
        r1 = True
    
    except:
        print("RF1 data error")
        r1 = False
    
    try:
        d2 = RFpower(path+"RFLog2.txt")
        d2.plotRF()
        r2 = True
    except:
        print("RF2 data error")
        r2 = False
    
    if r1 and r2:
        comboPlot(d1,d2)
    
    plt.show()
    
