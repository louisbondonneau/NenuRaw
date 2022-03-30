

from pathlib import Path
import numpy as np
import re
import os
import glob

from astropy.time import Time
import matplotlib
import matplotlib.pyplot as plt

from NenuRaw.luppi_structure import LUPPI
from NenuRaw.rawtf_structure import Rawtf
from NenuRaw.waveolaf_structure import WaveOlaf

# from multiprocessing import Pool, TimeoutError
import threading

from NenuRaw.wav_utils import Wav

# from pympler import summary
# from pympler import muppy
import time


class Raw:
    def __init__(self, obs_files, freq_start=0, freq_end=999, start=None, end=None, duration=None, block_start=0, block_end=1, verbose=False):
        self.freq_start, self.freq_end = freq_start, freq_end
        if(self.freq_end < self.freq_start):
            self.freq_end, self.freq_start = self.freq_start, self.freq_end

        self.block_start, self.block_end = block_start, block_end
        if(self.block_end < self.block_start):
            self.block_end, self.block_start = self.block_start, self.block_end

        self.names = obs_files
        self.fnames_raw = sorted(list(obs_files))
        self.data = []
        self.header = []
        self.verbose = verbose
        self.nof_blocks = []
        self.nof_subblocks = []
        self.fourier_methode = []
        print(self.fnames_raw)

        # LUPPI.print_header(LUPPI, fn_raw)
        # exit(0)
        fn_raw = Path(self.fnames_raw[0])
        if (LUPPI.try_luppi(LUPPI, fn_raw)):
            file_format = LUPPI(fn_raw, verbose=self.verbose)
        elif (WaveOlaf.try_waveolaf(fn_raw)):
            file_format = WaveOlaf(fn_raw, verbose=self.verbose)
        elif (Rawtf.try_rawtf(fn_raw)):
            file_format = Rawtf(fn_raw, verbose=self.verbose)
        else:
            raise NameError('444: No suitable format found')

        self.file_format = file_format

        if(start is not None):
            if (self.__try_iso(start)):
                start = self.__try_iso(start)
            elif (self.__try_isot(start)):
                start = self.__try_isot(start)
            else:
                start = float(start)
        else:
            start = 0.0

        if(end is not None):
            if (self.__try_iso(end)):
                end = self.__try_iso(end)
            elif (self.__try_isot(end)):
                end = self.__try_isot(end)
            else:
                end = float(end)
        elif(duration is not None):
            end = start + duration
        else:
            end = 1.0

        if(start > end):
            end, start = start, end

        self.requested_start = np.copy(start)
        self.requested_end = np.copy(end)

        for lane_idx, (fn_raw) in enumerate(self.fnames_raw):
            # print(lane_idx)
            print(fn_raw)
            fn_raw = Path(fn_raw)
            if verbose:
                print("Opening %s" % (fn_raw.name))

            # dt_header, header, dt_block = self.infer_data_structure(fn_raw)
            data = file_format.data

            self.data.append(data)  # keep it if multiple file in future
            self.header.append(file_format.header)

            # header = self.header[lane_idx]

            if (lane_idx == 0):  # if first file
                self.format_func = file_format
                self.nchan_file = file_format.nchan_file
                self.dm = file_format.dm
                self.nschan = file_format.nschan
                self.nschan_file = file_format.nschan_file
                self.nof_polcpx = file_format.nof_polcpx
                self.overlap = file_format.overlap
                self.tbin = file_format.tbin
                self.time_int = file_format.time_int
                self.beamletmin_file = file_format.beamletmin_file
                self.beamletmax_file = file_format.beamletmax_file
                self.chan_bw = file_format.chan_bw
                self.imjd = file_format.imjd
                self.smjd = file_format.smjd
                self.offmjd = file_format.offmjd

                self.time_off = 0

                if start is not None:
                    self.block_start = int(np.floor(start / self.time_int))
                    self.time_start = start
                else:
                    self.time_start = self.block_start * self.time_int
                    start = self.time_start

                if end is not None:
                    self.block_end = int(np.ceil(end / self.time_int))
                    self.time_end = end
                else:
                    self.time_end = self.block_end * self.time_int
                    end = self.time_end

                # subblock is needed?
                print("Block start %d Block end %d" % (self.block_start, self.block_end))
                self.subblocking = int(2**np.ceil(np.log(self.time_int / (self.time_end - self.time_start)) / np.log(2)))

                # print(self.time_end-self.time_start)
                # print(self.time_int/(self.time_end-self.time_start))
                # print(self.subblocking)

                if(self.subblocking < 2):
                    self.subblocking = 2
                self.block_start = int(np.floor(start / (self.time_int)))
                self.block_end = int(np.floor(end / (self.time_int)))

                self.subtime_int = self.time_int / float(self.subblocking)
                # self.subblock_start = int(np.floor(start/self.subtime_int)) % self.subblocking
                self.subblock_start = int(np.floor(start / self.subtime_int) - self.subblocking * self.block_start) % self.subblocking
                # self.subblock_end = int(np.ceil(end/self.subtime_int)) % self.subblocking
                self.subblock_end = int(np.floor(end / self.subtime_int) - self.subblocking * (self.block_end - 1)) % self.subblocking
                # if(self.subblock_end == 0):
                #    self.subblock_end = 1
                if(self.subblocking == 1):
                    self.subblock_end = int(1)
                # self.time_off += self.subblocking*self.time_int*np.floor(int(np.ceil(end/(self.time_int)))/self.subblocking)
                # self.time_off += self.subblock_start*self.subtime_int
                # self.nschan = int(self.nschan*(self.block_end-self.block_start)/self.subblocking)
                print("Subbloking by factor %d" % (self.subblocking))
                print("Block start %d Block end %d" % (self.block_start, self.block_end))
                print("subBlock start %d subBlock end %d" % (self.subblock_start, self.subblock_end))
                self.time_start = self.block_start * self.time_int + self.subblock_start * self.subtime_int
                self.time_end = (self.block_end) * self.time_int + (self.subblock_end + 1) * self.subtime_int

                print("Time stat %f. Time end %f" % (self.time_start, self.time_end))

                self.chan_bw_file = np.copy(self.chan_bw)
                # print(self.beamletmin, self.beamletmax, self.nchan)
                # print(lowchan, hichan, self.nchan)
                self.beamletmin, self.beamletmax = int(self.freq_start / self.chan_bw), int(self.freq_end / self.chan_bw)
                if (self.beamletmin < self.beamletmin_file) or (self.beamletmin > self.beamletmax_file):
                    self.beamletmin = self.beamletmin_file
                if (self.beamletmax > self.beamletmax_file) or (self.beamletmax < self.beamletmin_file):
                    self.beamletmax = self.beamletmax_file

                self.nchan = self.beamletmax_file - self.beamletmin_file + 1
                self.chan_factor = int(1)
                self.nsblk = self.nschan * self.nchan
                self.beamlet_k = np.linspace(self.beamletmin, self.beamletmax, self.nchan)
                self.minfreq = (self.beamletmin_file) * self.chan_bw
                self.maxfreq = (self.beamletmax_file) * self.chan_bw
                self.freq = np.linspace(self.minfreq, self.maxfreq, self.nchan)
                self.freq_file = np.copy(self.freq)
                self.mjd0 = self.imjd + (self.smjd + self.offmjd) / 86400.
            self.nof_blocks.append(file_format.nof_blocks)
            self.nof_subblocks.append([])
            # print((np.sum(self.nof_blocks)),  self.block_end)
            if ((np.sum(self.nof_blocks)) > self.block_end):
                self.nof_blocks[lane_idx] = int(self.block_end + 1 - np.sum(self.nof_blocks[:-1]))
            # print((np.sum(self.nof_blocks)),  self.block_end)

    def __try_isot(self, time_string):
        try:
            time_obj = Time(time_string, format='isot', scale='utc')  # '2013-12-04T00:00:00.666668'
        except:
            return False
        imjd = self.file_format.imjd
        smjd = self.file_format.smjd
        offmjd = self.file_format.offmjd
        mjd0 = imjd + (smjd + offmjd) / 86400.
        time_from_start = time_obj.mjd - mjd0
        if (time_from_start < 0):
            raise ('ERROR: input time %s is befor the start time %s' % (time_string, Time(mjd0, format='mjd', scale='utc').iso))
        # conv in sec from start
        sec_from_start = time_from_start * 86400.
        return sec_from_start

    def __try_iso(self, time_string):
        try:
            time_obj = Time(time_string, format='iso', scale='utc')  # '2013-12-04T00:00:00.666668'
        except:
            return False
        imjd = self.file_format.imjd
        smjd = self.file_format.smjd
        offmjd = self.file_format.offmjd
        mjd0 = imjd + (smjd + offmjd) / 86400.
        time_from_start = time_obj.mjd - mjd0
        if (time_from_start < 0):
            raise ('ERROR: input time %s is befor the start time %s' % (time_string, Time(mjd0, format='mjd', scale='utc').iso))
        # conv in sec from start
        sec_from_start = time_from_start * 86400.
        return sec_from_start

    def __nfft(self):
        # print("__nfft FFTLEN = %f" % self.fftlen)
        # print("__nfft NSCHAN = %f" % self.nschan)
        self.nfft = ((self.nschan / self.nof_polcpx)) / self.fftlen
        if (self.nfft <= 1):
            self.nfft = 1
            self.fftlen = int((self.nschan / self.nof_polcpx))
        self.nfft = int(self.nfft)

    def __ds(self):
        self.ds = self.ds_ms * 1e-3 / (self.tbin * self.fftlen)
        self.ds = int(2**np.floor(np.log(self.ds) / np.log(2)))
        if(self.ds == 0):
            self.ds = 1
        if(self.ds > self.nfft):
            self.ds = self.nfft
        if(self.nschan / self.nof_polcpx / self.fftlen % self.ds != 0):
            while (self.nschan / self.nof_polcpx / self.fftlen % self.ds != 0):
                self.ds -= 1
        self.ds_ms = float(self.ds) * self.tbin * float(self.fftlen) * 1e3

    def __df(self):
        self.df = int(2**round(np.log(self.df) / np.log(2)))
        self.df_mhz = float(self.df) * self.chan_bw / self.fftlen

    def fourier_computation(self, fftlen=16, ds_ms=30, df=1, file_num=0, pol="I"):
        self.ds_ms = ds_ms  # ms
        self.df = int(df)  # nchan

        # self.chan_sel_start = int(np.floor((self.beamletmin - self.beamletmin_file)/self.df))
        # self.chan_sel_stop = int(np.ceil(self.nchan_file - (self.beamletmax_file - self.beamletmax)))
        self.chan_sel_start = self.beamletmin - self.beamletmin_file
        self.chan_sel_stop = self.nchan_file - (self.beamletmax_file - self.beamletmax)

        self.fftlen = fftlen

        if (self.fftlen > (self.time_end - self.time_start) / self.tbin):
            self.fftlen = np.floor((self.time_end - self.time_start) / self.tbin)
            self.fftlen = int(2**np.floor(np.log(self.fftlen) / np.log(2))) / 2
            print("WARNING fftlen > nb sample per chan will use %d" % int(self.fftlen))

        self.fftlen = int(2**round(np.log(self.fftlen) / np.log(2)))
        self.__nfft()
        self.__ds()
        self.__df()

        print("  ds         is: %s fftlen %.3f ms" % (self.ds, self.ds_ms))
        print("  df         is: %s chan   %.3f MHz" % (self.df, self.df_mhz))
        print("  nfft       is: %s" % self.nfft)
        print("  fftlen     is: %s int %.3f ms" % (self.fftlen, self.fftlen * self.tbin * 1e3))
        # print(self.nsblk*5.12e-6/192/4)

        # ------ without mp ------
        for block_num in range(self.block_start, self.block_end + 1):
            data_plot = self.fourier_computation_block(block_num, file_num=file_num, pol=pol)
            if (block_num == int(self.block_start)):
                self.dynspec = data_plot
            else:
                self.dynspec = np.append(self.dynspec, data_plot, axis=2)

        # ------ try mp ------
        # pool = Pool(processes=int(8))
        # mp_processe_all = []
        # for block_id in range(self.block_start, self.block_end + 1):
        #     mp_processe = pool.apply_async(self.fourier_computation_block,
        #                                    (block_id, file_num, pol))
        #     mp_processe_all.append(mp_processe)
        # for i in mp_processe_all:
        #     i.get()
        # for block_id in range(len(mp_processe_all)):
        #     if (block_id == int(0)):
        #         self.dynspec = np.copy(mp_processe_all[block_id])
        #     else:
        #         self.dynspec = np.append(self.dynspec, mp_processe_all[block_id], axis=2)

        self.dynspec = np.flip(self.dynspec, axis=1)
        # print(np.sum(self.nof_subblocks))
        # print(np.shape(self.dynspec))
        # self.dynspec = np.reshape(self.dynspec, (int(self.nchan*self.fftlen/self.df), int((np.sum(self.nof_blocks)-self.block_start)*self.nfft/self.ds)))
        # print(np.sum(self.nof_blocks), self.block_start, self.block_end)
        # print(self.nschan_file, np.sum(self.nof_subblocks), self.nof_polcpx, self.subblocking, self.fftlen, self.ds)
        self.nschan = int(float(np.sum(self.nof_subblocks) * self.nschan_file) / self.nof_polcpx / self.subblocking / self.fftlen / self.ds)
        # print('ICCI', self.nschan)
        # nfft = self.nschan/self.nof_po
        self.dynspec = np.reshape(self.dynspec, (int((self.nchan / self.chan_factor) * self.fftlen / self.df), self.nschan))

        if(self.chan_sel_start != 0) or (self.chan_sel_stop != self.nchan - 1):
            self.nchan_file = np.copy(self.nchan)
            self.chan_sel_start, self.chan_sel_stop = 0, self.nchan - 1

        self.tbin /= self.fftlen * self.ds
        self.nof_blocks = int((np.sum(self.nof_blocks) - self.block_start))

        # if (self.subblocking > 1):
        #    self.nschan =  int(self.nschan/(self.fftlen*self.ds*self.nof_polcpx))
        # else:
        #    self.nschan =  self.nof_blocks * int(self.nschan/(self.fftlen*self.ds*self.nof_polcpx))
        self.nsblk = self.nschan * self.nchan / self.nof_blocks
        self.nof_polcpx = 1
        # self.ds = 1
        self.nchan = self.nchan * self.fftlen / self.df
        self.chan_factor = int(self.nchan / self.nchan_file)
        self.minfreq = (self.beamletmin - 0.5) * self.chan_bw_file + self.chan_bw / self.fftlen / 2
        self.maxfreq = (self.beamletmax + 0.5) * self.chan_bw_file - self.chan_bw / self.fftlen / 2
        self.chan_bw = self.chan_bw / self.fftlen / self.df
        self.freq = np.linspace(self.minfreq, self.maxfreq, int(self.fftlen * (self.beamletmax - self.beamletmin + 1) / self.df))

    def fourier_computation_block(self, block_num, file_num=0, pol='I', dm=None):
        if dm is not None:
            self.dm = dm
        if(self.verbose):
            print("dedisperse: dm = %.6f pc cm-3" % self.dm)
        # dt = self.dm * 4150 * (self.minfreq**-2 - self.maxfreq**-2.)
        # dt_bin = np.ceil(dt / self.tbin)

        # block_num  my_spectra.data mjd0 dat nchan nfft overlap fftlen nof_polcpx
        self.patidx = self.format_func.get_patidx(data=self.data, file_num=file_num, block_num=block_num)
        self.mjd = self.mjd0 + (self.patidx * self.tbin) / 86400.
        self.dates = Time(self.mjd, format='mjd', scale='utc').iso
        print('block: %d/%d start: %s pol: %s' % (block_num, (self.block_end), self.dates, pol))
        if (True):  # (self.subblocking > 1):
            if(block_num == self.block_start):
                self.dates_start = Time(self.mjd + float(self.subblock_start) * self.subtime_int / 86400., format='mjd', scale='utc').iso
                if(self.verbose):
                    print("block %d from subblock %d start: %s pol: %s" % (block_num, self.subblock_start, self.dates_start, pol))
            if(block_num == self.block_end):
                self.dates_end = Time(self.mjd + float(self.subblock_end + 1) * self.subtime_int / 86400., format='mjd', scale='utc').iso
                if(self.verbose):
                    print("block %d to subblock %d end: %s pol: %s" % (block_num, self.subblock_end, self.dates_end, pol))

        if (self.dm > 0):
            self.wav = self.__get_dedispers_block(block_num, file_num)
        else:
            self.wav = self.format_func.get_block(self.data, file_num, block_num)
            # self.wav = np.reshape(self.wav, (self.nchan_file, (int(self.nfft)+int(self.overlap/self.fftlen))*self.fftlen, self.nof_polcpx))
            self.wav = np.reshape(self.wav, (self.nchan_file, int(self.nschan_file / self.nof_polcpx + self.overlap), self.nof_polcpx))
            if (self.overlap > 0):
                self.wav = self.wav[:, :-self.overlap, :]
        # print(type(self.data[file_num]['data'][0,0]))
        # print(type(self.data[file_num]['data'][block_num][0]));
        # print(type(self.wav[0, 0, 0])); exit(0)

        # self.wav = self.wav.astype('float32').view('complex64')

        self.wav.shape = (self.nchan_file, self.subblocking, int(self.nschan_file / self.nof_polcpx / self.subblocking), self.nof_polcpx)
        self.wav = np.reshape(self.wav, (self.nchan_file, self.subblocking, int(self.nschan_file / self.nof_polcpx / self.subblocking), self.nof_polcpx))
        # print(np.shape(self.wav))
        if (self.block_start == self.block_end):
            self.wav = self.wav[:, self.subblock_start:self.subblock_end + 1, :, :]
            self.nschan = int(self.nschan_file * (self.subblock_end + 1 - self.subblock_start) / self.subblocking)
            self.nof_subblocks[file_num].append(int(self.subblock_end + 1 - self.subblock_start))
        elif(block_num == self.block_start):
            self.wav = self.wav[:, self.subblock_start:, :, :]
            self.nschan = int(self.nschan_file * (self.subblocking - self.subblock_start) / self.subblocking)
            self.nof_subblocks[file_num].append(int(self.subblocking - self.subblock_start))
        elif(block_num == self.block_end):
            self.wav = self.wav[:, :self.subblock_end + 1, :, :]
            self.nschan = int(self.nschan_file * (self.subblock_end + 1) / self.subblocking)
            self.nof_subblocks[file_num].append(int(self.subblock_end + 1))
        else:
            self.nschan = int(self.nschan_file)
            self.nof_subblocks[file_num].append(int(self.subblocking))
        self.blocksize = self.nschan * self.nchan_file
        # print('LAA: ', self.nchan_file, self.nschan, self.nof_polcpx, self.nof_polcpx)
        self.wav = np.reshape(self.wav, (self.nchan_file, int(self.nschan / self.nof_polcpx), self.nof_polcpx))

        # if( self.overlap > 0 ):
        #    self.wav = self.wav[:, :int(self.nschan_file/self.nof_polcpx - self.overlap), :]
        #    self.nschan_file = int(self.nschan_file - self.overlap*self.nof_polcpx)
        #    self.nschan = self.nschan_file
        #    self.blocksize = self.nschan * self.nchan_file
        #    self.overlap = 0

        self.__nfft()
        # self.__ds()

        self.wav = self.wav[self.chan_sel_start:self.chan_sel_stop, :, :]
        self.nchan = self.chan_sel_stop - self.chan_sel_start

        self.wav = self.wav.astype('float32').view('complex64')
        npol = int(self.nof_polcpx / 2)

        for fourier_func in self.fourier_methode:
            print(fourier_func)
            fourier_func(self)

        self.wav.shape = (self.nchan, int(self.nfft), self.fftlen, npol)
        self.wav = np.reshape(self.wav, (self.nchan, int(self.nfft), self.fftlen, npol))

        # print('0:',np.shape(self.wav))
        self.wav = self.wav[:, :int(self.nfft), :, :]
        # self.wav = self.wav.astype('float32').view('complex64')
        # print('1:',np.shape(self.wav))
        if(self.fftlen > 1):
            self.wav = np.fft.fft(self.wav, axis=2)
            self.wav = np.fft.fftshift(self.wav, axes=2)  # ifftshift
        # print('2:',np.shape(self.wav))
        self.wav_X = self.wav[:, :, :, 0]
        self.wav_Y = self.wav[:, :, :, 1]

        # Xconj = self.wav_X.real**2 + self.wav_X.imag**2
        # Yconj = self.wav_Y.real**2 + self.wav_Y.imag**2
        # XYconj = self.wav_X*(self.wav_Y.real - self.wav_Y.imag)
        # YXconj = self.wav_Y*(self.wav_X.real - self.wav_X.imag)
        # I = Xconj + Yconj
        # Q =  Xconj - Yconj
        # U = 2*XconjY.real
        # V = 2*XconjY.imag
        # data_plot = self.wav_X.real**2 + self.wav_X.imag**2 + self.wav_Y.real**2 + self.wav_Y.imag**2
        if (pol == "Q"):
            Xconj = self.wav_X.real**2 + self.wav_X.imag**2
            Yconj = self.wav_Y.real**2 + self.wav_Y.imag**2
            data_plot = Xconj - Yconj
        elif (pol == "U"):
            XconjY = self.wav_X.conj() * self.wav_Y
            data_plot = 2 * XconjY.real
        elif (pol == "L"):
            Xconj = self.wav_X.real**2 + self.wav_X.imag**2
            Yconj = self.wav_Y.real**2 + self.wav_Y.imag**2
            XconjY = self.wav_X.conj() * self.wav_Y
            data_plot = (Xconj - Yconj)**2 + (2 * XconjY.real)**2
        elif (pol == "V"):
            XconjY = self.wav_X.conj() * self.wav_Y
            data_plot = 2 * XconjY.imag
        elif (pol == "XX"):
            XconjX = self.wav_X.conj() * self.wav_X
            data_plot = XconjX
        elif (pol == "YY"):
            YconjY = self.wav_Y.conj() * self.wav_Y
            data_plot = YconjY
        else:  # "I"
            Xconj = self.wav_X.real**2 + self.wav_X.imag**2
            Yconj = self.wav_Y.real**2 + self.wav_Y.imag**2
            data_plot = Xconj + Yconj

        # print('3:',np.shape(data_plot))
        # data_plot = np.reshape(data_plot, (self.nchan, int(self.nfft/self.ds), self.ds, self.fftlen))
        # print(self.nfft, self.ds, float(self.nfft/self.ds))
        # print(self.nchan,  int(self.nfft/self.ds), self.ds, int(self.fftlen/self.df), self.df)
        data_plot = np.reshape(data_plot, (self.nchan, int(self.nfft / self.ds), self.ds, int(self.fftlen / self.df), self.df))

        # data_plot = data_plot[ self.chan_sel_start: self.chan_sel_stop]

        data_plot = np.mean(data_plot, axis=4)
        data_plot = np.mean(data_plot, axis=2)

        # data_plot = np.reshape(data_plot, (nchan, ds, int(nfft/ds), fftlen))
        # data_plot = np.mean(data_plot, axis=1)
        data_plot = np.rot90(data_plot, k=1, axes=(1, 2))
        # np.rot90(data_plot, k=1, axes=(1, 2))
        # print('4:',np.shape(data_plot))
        return data_plot

        # allObjects = muppy.get_objects()
        # sum = summary.summarize(allObjects)
        # summary.print_(sum)
        # time.sleep(10)
        # exit(0)

    def __get_dedispers_block(self, block_num, file_num):
        tmp = np.copy(self.format_func.get_block(self.data, file_num, block_num))

        blocksize = float(self.nchan_file * self.nschan_file / self.nof_polcpx)

        tmp.shape = (self.nchan_file, int(blocksize / self.nchan_file + self.overlap), self.nof_polcpx)
        tmp = np.reshape(tmp, (self.nchan_file, int(blocksize / self.nchan_file + self.overlap), self.nof_polcpx))

        if(self.overlap > 0):
            tmp = tmp[:, :-self.overlap, :]

        blocksize_perchan = float(self.nschan_file / self.nof_polcpx)
        old_block_num = block_num - 1

        for ichan in range(self.nchan_file - 2, -1, -1):
            dt = self.dm * 4150. * (self.freq[ichan]**-2 - self.freq[-1]**-2.)
            start_bin = np.round(dt / self.tbin)
            end_bin = start_bin + blocksize_perchan - 1
            start_block = int(np.floor(start_bin / blocksize_perchan))
            end_block = int(np.floor(end_bin / blocksize_perchan))
            start_bin = int(start_bin - (start_block * blocksize_perchan))
            end_bin = int(end_bin - (end_block * blocksize_perchan))

            if (block_num + start_block != old_block_num):
                current_block = self.format_func.get_block(self.data, file_num, block_num + start_block)
                current_block.shape = (self.nchan_file, int(blocksize / self.nchan_file + self.overlap), self.nof_polcpx)
                current_block = np.reshape(current_block, (self.nchan_file, int(blocksize / self.nchan_file + self.overlap), self.nof_polcpx))
                if(self.overlap > 0):
                    current_block = current_block[:, :-self.overlap, :]
                current_block_next = self.format_func.get_block(self.data, file_num, block_num + end_block)
                current_block_next.shape = (self.nchan_file, int(blocksize / self.nchan_file + self.overlap), self.nof_polcpx)
                current_block_next = np.reshape(current_block_next, (self.nchan_file, int(blocksize / self.nchan_file + self.overlap), self.nof_polcpx))
                if(self.overlap > 0):
                    current_block_next = current_block_next[:, :-self.overlap, :]

                old_block_num = block_num + start_block
            # print(np.shape(tmp))
            # print(start_bin, end_bin)
            tmp[ichan, :int(blocksize_perchan - start_bin), :] = np.copy(current_block[ichan, start_bin:, :])
            if(end_bin != blocksize_perchan - 1):
                tmp[ichan, int(blocksize_perchan - start_bin + 1):, :] = np.copy(current_block_next[ichan, :end_bin, :])

            # print(self.freq[ichan], start_block, start_bin, end_block, end_bin)
        # tmp.shape = (self.nchan_file, int(self.nfft)+int(self.overlap/self.fftlen), self.fftlen, self.nof_polcpx)
        # tmp = np.reshape(tmp, (self.nchan_file, int(self.nfft)+int(self.overlap/self.fftlen), self.fftlen, self.nof_polcpx))
        return tmp

    def new_fourier_methode(self, func):
        self.fourier_methode.append(func)

    def print_fourier_methode(self):
        print("")
        for i in range(len(self.fourier_methode)):
            print("methode %d:  %s" % (i, self.fourier_methode[i]))


class Thread (threading.Thread):
    def __init__(self, func, args=((),), kwargs={}, name=None):
        threading.Thread.__init__(self)  # init mother class
        self.func2thread = func
        self.args = args
        self.kwargs = kwargs
        self.name = name

    def run(self):
        self.func2thread(*self.args, **self.kwargs)
