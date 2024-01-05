import numpy as np
import matplotlib.pyplot as plt
import sys

class OceanSpectra():

    def __init__(self,fin):

        self.fname = fin
        self.loadData(fin)


        # for plotting peaks
        self.lines = []
        self.freqs = []


    def loadData(self,fin):

        with open(fin) as f:
            indata = f.readlines()


        # get line where data begins
        j_start = np.argwhere([line.find(">>>")+1 for line in indata])[0][0]

        # parse everything prior into meta data dict
        meta = {} 
        for line in indata[2:j_start]:

            string = line.strip()
            k = string.find(':')

            key = string[:k]
            val = string[k+2:]
            meta[key] = val

        # parse spectral axis (nm)
        spectral_axis = np.array(indata[j_start+1].strip().split('\t'),float)

        # parse spectra
        human_time = []
        unix_time = []
        data = []
        for line in indata[j_start+2:]:
            stream = line.strip().split('\t')

            human_time.append(stream[0])
            unix_time.append(stream[1])
            data.append(stream[2:])


        unix_time = np.array(unix_time, int)
        data = np.array(data, float)


        # save
        self.raw_data = indata
        self.human_time = human_time
        self.unix_time = unix_time
        self.data = data
        self.meta = meta

        self.N_spectra = len(data)

        # get time
        self.dt = float(meta['Integration Time (sec)'])
        self.time_axis = unix_time - unix_time[0] # ms
        self.spectral_axis = spectral_axis

    def findPeak(self,f0=656.363, # nm
                 ):

        # H-alpha 656.363 approx
        # H-beta 486.0675 approx

        freq = self.spectral_axis
        j0 = np.argmin( np.abs(freq - f0) )

        f_time = self.data[:,j0] 

        self.lines.append(f_time)
        self.freqs.append(f0)

    def plotPeaks(self, axs):
        # plot identified freq peaks over time
        time = self.time_axis

        N_lines = len(self.lines)
        for j in np.arange(N_lines):
            axs.plot(time, self.lines[j], label=f"{self.freqs[j]} nm")

    def plot2d(self, j=50, save=False):

        s_ax = self.spectral_axis
        t_ax = self.time_axis
        data = self.data

        fig,axs = plt.subplots(2,1,figsize=(10,5))

        C = axs[0].contourf(t_ax, s_ax, data.T, cmap='inferno')
        axs[0].set_ylabel('wavelength (nm)')
        axs[0].set_xlabel('time (ms)')

        t_slice = t_ax[j]
        axs[0].axvline(t_slice, color='r', ls='--')
        fig.colorbar(C)

        axs[1].plot( s_ax, data[j], 'r', lw=0.7, label=f"t = {t_slice} ms" )
        axs[1].set_xlabel('wavelength (nm)')
        axs[1].set_ylabel('counts')
        axs[1].legend()
        axs[1].grid()

        fig.suptitle(self.fname)
        fig.tight_layout()

        if save:
            fig.savefig(save)


### Test Driver
if __name__ == "__main__":
    fin = sys.argv[1]
    data = OceanSpectra(fin)

    data.plot2d()
    plt.show()
