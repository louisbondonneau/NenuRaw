

from pathlib import Path
import numpy as np
import re
import os

from astropy.time import Time
import matplotlib
import matplotlib.pyplot as plt

from NenuRaw.raw_utils import Raw

#from pympler import summary
#from pympler import muppy
#import time


class Dynspec(Raw):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def dedisperse(self, dm=None):
        if dm is not None:
            self.dm = dm
            print(r"dedisperse: DM = %.6f pc cm$^{-3}$" % self.dm)
            if (self.dm > 0):
                dt = self.dm * 4150 * (self.minfreq**-2 - self.maxfreq**-2.)
                dt_bin = np.round(dt / self.tbin)  # ceil -> round
                self.time_off = -int(dt_bin) * self.tbin
                pad = np.nan * np.zeros([int(self.nchan), int(dt_bin)])
                self.dynspec = np.concatenate((pad, self.dynspec), axis=1)
                for ichan in range(len(self.freq)):
                    dt = self.dm * 4150 * (self.freq[ichan]**-2 - self.maxfreq**-2.)
                    dt_bin = np.ceil(dt / self.tbin)
                    self.dynspec[ichan, :] = np.roll(self.dynspec[ichan, :], -int(dt_bin))
            elif(self.dm < 0):
                print("DM must be possitive")

    def clean(self, threshold=3):
        median = np.median(self.dynspec)
        MAD = np.median(np.abs(self.dynspec - median))
        ind0 = np.where(self.dynspec > (median + threshold * MAD))
        ind1 = np.where(self.dynspec < (median - threshold * MAD))
        self.dynspec[ind0] = median + threshold * MAD
        self.dynspec[ind1] = median - threshold * MAD

    def rm_dm0(self):
        time_serie = np.median(self.dynspec, axis=0)
        median = np.median(time_serie)
        for ichan in range(self.nchan):
            self.dynspec[ichan, :] += median - time_serie

    def rm_baseline(self):
        nchan, nsub = np.shape(self.dynspec)
        spectrum = np.mean(self.dynspec, axis=1)
        spectrum_file = np.mean(np.reshape(spectrum, (self.nchan_file, int(self.chan_factor))), axis=1)
        median = np.median(self.dynspec)
        for ichan in range(self.nchan_file):
            # for ifft in range(int((np.sum(self.nof_blocks)-self.block_start)*self.nfft/self.ds)):
            # for ifft in range(nsub-1):
            for ifft in range(self.nschan):
                self.dynspec[int(ichan * self.chan_factor):int((ichan + 1) * self.chan_factor), ifft] = self.dynspec[int(ichan * self.chan_factor)                                                                                                                     :int((ichan + 1) * self.chan_factor), ifft] - spectrum[int(ichan * self.chan_factor):int((ichan + 1) * self.chan_factor)]
                self.dynspec[int(ichan * self.chan_factor):int((ichan + 1) * self.chan_factor), ifft] = self.dynspec[int(ichan *
                                                                                                                         self.chan_factor):int((ichan + 1) * self.chan_factor), ifft] * (median / spectrum_file[ichan]) + median

    def rm_pfb(self):
        # (192*8, 1024)
        # self.dynspec = np.reshape(self.dynspec, (self.nchan_file, self.chan_factor, int((np.sum(self.nof_blocks)-self.block_start)*self.nfft/self.ds)))
        # (192, 8, 1024)
        spectrum = np.mean(self.dynspec, axis=1)
        self.pfb = np.median(np.reshape(spectrum, (self.nchan_file, self.chan_factor)), axis=0)
        self.pfb /= np.median(self.pfb)
        # PFB = np.median(PFB, axis=0)
        # median = np.median(self.dynspec)
        for ichan in range(self.nchan_file):
            for ifft in range(self.nschan):
                median = np.median(self.dynspec[ichan * self.chan_factor:(ichan + 1) * self.chan_factor, ifft])
                self.dynspec[ichan * self.chan_factor:(ichan + 1) * self.chan_factor, ifft] = self.dynspec[ichan * self.chan_factor:(
                    ichan + 1) * self.chan_factor, ifft] - (self.pfb - 1) * median
        # self.dynspec = np.reshape(self.dynspec, (self.nchan_file*self.chan_factor, int((np.sum(self.nof_blocks)-self.block_start)*self.nfft/self.ds)))

    def plot_dynspec(self):
        if ((self.time_end - self.time_start) > 60):
            time_ratio = 1 / 60.
            label_start = (self.time_start + self.time_off) * time_ratio
            diff = np.floor(label_start)
            ref_time = self.mjd0 + (diff / time_ratio / 86400.)
            label_start = label_start - diff
            label_end = (self.time_end + self.time_off) * time_ratio - diff
            xlabel = 'Time (min) from %s' % (Time(ref_time, format='mjd', scale='utc').iso)
        elif ((self.time_end - self.time_start) > 1):
            time_ratio = 1.
            label_start = (self.time_start + self.time_off) * time_ratio
            diff = np.floor(label_start)
            ref_time = self.mjd0 + (diff / time_ratio / 86400.)
            label_start = label_start - diff
            label_end = (self.time_end + self.time_off) * time_ratio - diff
            xlabel = 'Time (sec) from %s' % (Time(ref_time, format='mjd', scale='utc').iso)
        elif ((self.time_end - self.time_start) > 0.001):
            time_ratio = 1.0e3
            label_start = (self.time_start + self.time_off) * time_ratio
            diff = np.floor(label_start)
            ref_time = self.mjd0 + (diff / time_ratio / 86400.)
            label_start = label_start - diff
            label_end = (self.time_end + self.time_off) * time_ratio - diff
            xlabel = 'Time (ms) from %s' % (Time(ref_time, format='mjd', scale='utc').iso)
        else:
            time_ratio = 1.0e6
            label_start = (self.time_start + self.time_off) * time_ratio
            diff = np.floor(label_start)
            ref_time = self.mjd0 + (diff / time_ratio / 86400.)
            label_start = label_start - diff
            label_end = (self.time_end + self.time_off) * time_ratio - diff
            xlabel = 'Time (microsec) from %s' % (Time(ref_time, format='mjd', scale='utc').iso)
        plt.subplots_adjust(hspace=0, wspace=0)
        ax0 = plt.subplot2grid((4, 1), (0, 0), colspan=1, rowspan=1)
        time_serie = np.nanmean(self.dynspec, axis=0)
        time_vector = np.linspace(label_start, label_end, len(time_serie))
        ax0.plot(time_vector, time_serie)
        ax0.axes.get_xaxis().set_visible(False)
        ax0.set_title(os.path.basename(self.names[0]))
        ax1 = plt.subplot2grid((4, 1), (1, 0), colspan=1, rowspan=3, sharex=ax0)
        ax1.imshow(np.flipud(self.dynspec), interpolation='none', cmap='afmhot', aspect='auto', extent=[
                   label_start, label_end, self.minfreq - self.chan_bw / 2, self.maxfreq + self.chan_bw / 2])
        ax1.set_ylabel('Frequency (MHz)')

        ax1.set_xlabel(xlabel)
        ax1.set_xlim([self.requested_start * time_ratio - diff, self.requested_end * time_ratio - diff])
        plt.show()

    def plot_spectra(self):
        spectrum = np.nanmean(self.dynspec, axis=1)
        PFB = np.reshape(spectrum, (self.nchan_file, self.chan_factor))
        PFB = np.median(PFB, axis=0)
        PFB = np.tile(PFB, int(self.nchan_file)) - np.mean(PFB)

        spectrum_mean = np.median(spectrum)

        for ichan in range(int(self.nchan_file)):
            ratio = np.median(spectrum[ichan * self.chan_factor: (ichan + 1) * self.chan_factor]) / spectrum_mean
            PFB[ichan * self.chan_factor: (ichan + 1) * self.chan_factor] *= ratio

        ax0 = plt.subplot2grid((1, 1), (0, 0), colspan=1, rowspan=1)
        ax0.plot(self.freq, spectrum, label='data')
        ax0.plot(self.freq, spectrum - PFB, label='data-PFB')
        ax0.plot(self.freq, PFB, '--', label='PFB')
        ax0.legend(loc='upper right')
        ax0.grid(True, which="both", ls="-", alpha=0.65)
        ax0.set_xlabel('Frequency (MHz)')
        ax0.set_ylabel('Amplitude (AU)')
        ax0.set_title(os.path.basename(self.names[0]))
        plt.show()

    def print_info(self):
        start_string = Time(self.mjd0, format='mjd', scale='utc').isot
        print("\nstart time is  %s (isot)" % (start_string))
        print("max frequency  %f (MHz)" % (self.maxfreq))
        print("min frequency  %f (MHz)" % (self.minfreq))
        print("block duration %f (sec)\n" % (self.time_int))

    def get_spectra(self):
        print("WARNING: get_spectra() to the benefit of get_data()")
        # spectrum = np.nanmean(self.dynspec, axis=1)
        return self.dynspec

    def get_data(self):
        return self.dynspec

    def get_freqarray(self):
        return np.linspace(self.minfreq, self.maxfreq, int(self.nchan))

    def get_timearray(self):
        time_array = (np.linspace(self.time_start + self.time_off, self.time_end + self.time_off, int(self.nschan)) / 86400.) + self.mjd0
        time_array = Time(time_array, format='mjd', scale='utc')
        return time_array

    def get_start(self):
        return Time(self.mjd0, format='mjd', scale='utc').iso
