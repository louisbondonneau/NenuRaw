
import numpy as np


class Wav:
    def __init__(self):
        print("initialisation")

    def coherent_dedisp(self, raw):
        print("coherent_dedisp")
        # raw.dm = 2
        if(raw.dm != 0):
            dmfac = raw.dm * 2.0 * np.pi / 2.41e-10  # (2.41e-10*(1.0+s->earth_z4/1.0e4)) doppler TODO
            nchan, max_fftlen, npol = np.shape(raw.wav)
            self.__fft(raw, fftlen=max_fftlen)
            chirp = np.zeros(max_fftlen) + 1j * np.zeros(max_fftlen)
            for ichan in range(nchan):
                max_fftlen_vec = np.linspace(0, max_fftlen - 1, max_fftlen)
                # for i in range(max_fftlen):
                dfreq = max_fftlen_vec * raw.chan_bw_file / max_fftlen
                dfreq[(max_fftlen_vec > max_fftlen / 2 - 1)] -= raw.chan_bw_file  # if not use fftshift
                freqfac = dfreq / raw.freq[int(raw.chan_sel_start + ichan)]
                freqfac = freqfac * freqfac / (raw.freq[int(raw.chan_sel_start + ichan)] + dfreq)
                arg = dmfac * freqfac
                chirp[:].real = (np.cos(arg) / max_fftlen)
                chirp[:].imag = -1.0 * (np.sin(arg) / max_fftlen)
                for ipol in range(npol):
                    wav = np.copy(raw.wav[ichan, :, :, ipol]) * max_fftlen / raw.nschan_file / raw.nof_polcpx
                    raw.wav[ichan, :, :, ipol].real = wav.real * chirp.real - wav.imag * chirp.imag
                    raw.wav[ichan, :, :, ipol].imag = wav.imag * chirp.real + wav.real * chirp.imag
            self.__ifft(raw)
        else:
            print("WARNING: no need to coherently dedisperse at DM=0")

    def coherent_defaraday(self, raw):
        print("coherent_defaraday")

    def wav_cleaning_freq(self, raw):
        print("wav_cleaning_freq")
        nchan, max_fftlen, npol = np.shape(raw.wav)
        self.__fft(raw, fftlen=max_fftlen)
        self.__clean(raw)  # Hz cleaning
        self.__ifft(raw)

    def __clean(self, raw):
        sigma = 3
        nchan, nschan_perpol, fftlen, npol = np.shape(raw.wav)

        wav_X = raw.wav[:, :, :, 0]
        wav_Y = raw.wav[:, :, :, 1]

        Xconj = wav_X.real**2 + wav_X.imag**2
        Yconj = wav_Y.real**2 + wav_Y.imag**2
        full_intensity = Xconj + Yconj

        for ichan in range(nchan):
            for isub in range(nschan_perpol):
                current_fft = np.abs(full_intensity[ichan, isub, :])
                median = np.median(current_fft)
                std = 1.48 * (np.median(np.abs(current_fft - median)))
                bad_id = (current_fft > median + sigma * std)
                raw.wav[ichan, isub, bad_id, :] = 0 + 0j

    def wav_cleaning_time(self, raw):
        print("wav_cleaning_time")
        nchan, max_fftlen, npol = np.shape(raw.wav)
        print("clean_time")
        sigma = 3

        wav_X = raw.wav[:, :, 0]
        wav_Y = raw.wav[:, :, 1]

        Xconj = wav_X.real**2 + wav_X.imag**2
        Yconj = wav_Y.real**2 + wav_Y.imag**2
        full_intensity = Xconj + Yconj

        for ichan in range(nchan):
            current_fft = np.abs(full_intensity[ichan, :])
            median = np.median(current_fft)
            std = 1.48 * (np.median(np.abs(current_fft - median)))
            bad_id = (current_fft > median + sigma * std)
            raw.wav[ichan, bad_id, :] = 0 + 0j

    def __fft(self, raw, fftlen):
        # print("fft")
        nchan, nschan_perpol, npol = np.shape(raw.wav)
        nfft = nschan_perpol / fftlen
        # print('0:',np.shape(raw.wav))
        raw.wav.shape = (nchan, int(nfft), fftlen, npol)
        raw.wav = np.reshape(raw.wav, (nchan, int(nfft), fftlen, npol))
        # print('1:',np.shape(raw.wav))
        raw.wav = np.fft.fft(raw.wav, axis=2)
        # raw.wav = np.fft.fftshift(raw.wav, axes=2) #ifftshift
        # print('2:',np.shape(raw.wav))

    def __ifft(self, raw):
        # print("ifft")
        nchan, nschan_perpol, fftlen, npol = np.shape(raw.wav)
        # raw.wav = np.fft.ifftshift(raw.wav, axes=2)
        raw.wav = np.fft.ifft(raw.wav, axis=2)
        # raw.wav.shape = (nchan, nschan_perpol*fftlen, npol)
        raw.wav = np.reshape(raw.wav, (nchan, nschan_perpol * fftlen, npol))
        # print('3:',np.shape(raw.wav))
