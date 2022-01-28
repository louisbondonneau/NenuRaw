from NenuRaw import dynspec_utils
from NenuRaw import wav_utils

if __name__ == "__main__":
    # wavfile (GUPPI or RAWTF format)
    files = ['/data/lbondonneau/test/RAW/B1508+55_5586.rawolaf']

    # initialisation of the Dynspec object containing the methodes
    my_spectra = dynspec_utils.Dynspec(files,
                                       verbose=False,
                                       freq_start=0,
                                       freq_end=99,
                                       start=10,  # start time in sec
                                       end=10.4  # stop time in sec
                                       )

    # force dm to 0 pc cm-3 or will use the dm from the header in PSR obs
    my_spectra.dm = 0

    # execution of the processing method in the Fourier domain creating 64 sub-channels (~3 kHz channel)
    # and conversion to total intensity with time integration of 100 ms
    my_spectra.fourier_computation(fftlen=1, ds_ms=4 * 5.12e-6, pol="I")  # I, Q, U, L, V, XX, YY

    # plot spectrum
    # my_spectra.clean(threshold=10) # clean with threshold=10
    my_spectra.plot_spectra()
