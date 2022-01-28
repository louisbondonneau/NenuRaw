#!/usr/bin/env python3

from NenuRaw import dynspec_utils
from NenuRaw import wav_utils

if __name__ == "__main__":
    # wavfile (GUPPI or RAWTF format)
    files = ['/databf2/nenufar-pulsar/ES03/2020/05/B0834+06_D20200510T1700_58979_250507_0069_BEAM1.0000.raw']

    # initialisation of my_spectra object containing the methodes
    my_spectra = dynspec_utils.Dynspec(files,
                                       verbose=True,
                                       freq_start=0,  # min freq allow
                                       freq_end=99,  # max freq allow
                                       start=10,  # start time in sec
                                       end=20  # stop time in sec
                                       )

    # will use the dm from the header (can be changed with my_spectra.dm = X)

    # initialisation of my_wav_obj object containing the Fourier methodes as coherent dedispersion
    my_wav_obj = wav_utils.Wav()
    my_spectra.new_fourier_methode(my_wav_obj.coherent_dedisp)

    # execution of the processing method in the Fourier domain creating 1 sub-channels (196 kHz channel)
    # and conversion to total intensity with time integration of 5.12 microsec
    my_spectra.fourier_computation(fftlen=1, ds_ms=10, pol="I")  # I, Q, U, L, V, XX, YY

    my_spectra.clean(threshold=10)  # first clean with threshold=10
    my_spectra.rm_baseline()  # baseline flatening
    my_spectra.clean(threshold=4)  # second and refined cleaning
    my_spectra.rm_baseline()  # baseline refined flatening

    # plot of the dynspectrum
    my_spectra.plot_dynspec()

    # DIY methodes to play with the data
    my_dynspec_data = my_spectra.get_data()
    date = my_spectra.get_start()
    print(date)
    my_freq_array = my_spectra.get_freqarray()
    my_time_array = my_spectra.get_timearray()
    ###
