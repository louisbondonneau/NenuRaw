#!/usr/bin/env python3

from NenuRaw import dynspec_utils
from NenuRaw import wav_utils

if __name__ == "__main__":
    # wavfile (GUPPI or RAWTF format)
    files = ['/databf2/nenufar-pulsar/ES03/2020/05/B0834+06_D20200510T1700_58979_250507_0069_BEAM1.0000.raw']

    # initialisation of the Dynspec object containing the methodes
    my_spectra = dynspec_utils.Dynspec(files,
                                       verbose=True,
                                       freq_start=0,
                                       freq_end=99,
                                       start='2020-05-10T17:12:05.666668',  # start time in iso format or sec from start
                                       duration=2.5  # duration time in sec
                                       )

    # force dm to 0 pc cm-3 or will use the dm from the header in PSR obs
    my_spectra.dm = 0

    # execution of the processing method in the Fourier domain creating 64 sub-channels (~3 kHz channel)
    # and conversion to total intensity with time integration of 100 ms
    my_spectra.fourier_computation(
        fftlen=64, ds_ms=100, pol="I")  # I, Q, U, L, V, XX, YY

    # plot spectrum
    my_spectra.clean(threshold=10)  # clean with threshold=10
    my_spectra.plot_spectra()

    # DIY methodes to play with the data
    date = my_spectra.get_start()
    print(date)
    my_spectra.rm_pfb()  # clean spectra from the pfb
    # get_spectra will scrunch in time will get_data will keep the time domain
    my_dynspec_data = my_spectra.get_data()
    my_freq_array = my_spectra.get_freqarray()
    timearray = my_spectra.get_timearray()
    print(timearray.iso)
    ###
