#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
from NenuRaw_module import Dynspec
from NenuRaw_module import Wav

if __name__ == "__main__":
    # Setup argparse
    parser = argparse.ArgumentParser(description="Dynamic spectrum processing script.")
    parser.add_argument('--files', nargs='+', default=['/databf/nenufar-pulsar/DATA/ZENITH/RAW/ZENITH_D20231121T1436_60269_254191_0057_BEAM1.0000.raw'], help="List of files to process")
    parser.add_argument('--start', type=str, default="10", help="Start time in seconds from the start of the recording")
    parser.add_argument('--duration', type=float, default=20, help="Duration time in seconds")
    parser.add_argument('--fftlen', type=int, default=131072, help="FFT length for Fourier computation")
    parser.add_argument('--ds_ms', type=float, default=100.0, help="Time integration in milliseconds")
    parser.add_argument('--dm', type=float, default=0.0, help="Time integration in milliseconds")
    parser.add_argument('--pol', type=str, default='I', help="polarization in 'I' 'Q' 'U' 'L' 'V' 'XX' 'YY'")
    parser.add_argument('--freq_start', type=float, default=43.5, help="Start frequency in MHz")
    parser.add_argument('--freq_end', type=float, default=46.5, help="End frequency in MHz")
    parser.add_argument('--clean', action='store_true', help="End frequency in MHz")

    parser.add_argument('--clean_wave_freq', action='store_true', help="End frequency in MHz")
    parser.add_argument('--clean_wave_time', action='store_true', help="End frequency in MHz")

    parser.add_argument('--rm_baseline', action='store_true', help="End frequency in MHz")

    parser.add_argument('--plot_spectrum', action='store_true', help="End frequency in MHz")
    parser.add_argument('--plot_dynspec', action='store_true', help="End frequency in MHz")

    # Parse arguments
    args = parser.parse_args()
    num_files = len(args.files)
    # Créer une grille de subplots
    fig, axs = plt.subplots(num_files, 1, figsize=(10, num_files * 5), sharex=True, sharey=True)


    # Process each file in the list
    for index, file in enumerate(args.files):
        my_spectra = Dynspec([file],
                             verbose=True,
                             freq_start=args.freq_start,
                             freq_end=args.freq_end,
                             start=args.start,
                             duration=args.duration,
                             )
        
        my_spectra.dm = args.dm
        if (args.dm != 0) or (args.clean_wave_freq) or (args.clean_wave_time):
            my_wav_obj = Wav()
            if (args.dm != 0):
                my_spectra.new_fourier_methode(my_wav_obj.coherent_dedisp)
            if (args.clean_wave_freq):
                my_spectra.new_fourier_methode(my_wav_obj.wav_cleaning_freq)
            if (args.clean_wave_time):
                my_spectra.new_fourier_methode(my_wav_obj.wav_cleaning_time)

        my_spectra.fourier_computation(fftlen=args.fftlen,
                                       ds_ms=args.ds_ms,
                                       pol=args.pol
                                       )
        
        if (args.clean):
            my_spectra.clean(threshold=30)
        if (args.rm_baseline):
            my_spectra.rm_baseline()
        if (args.clean):
            my_spectra.clean(threshold=10)
            if (args.rm_baseline):
                my_spectra.rm_baseline()
        if (args.plot_dynspec):
            my_spectra.rm_pfb()
            if len(args.files) == 1:
                my_spectra.plot_dynspec(log10=True)
            else:
                my_spectra.plot_dynspec()
        
        if (args.plot_spectrum):
            if len(args.files) == 1:
                #my_spectra.plot_spectra() # only if single file
                my_spectra.rm_pfb()
                my_dynspec_data = my_spectra.get_data() # (nsubint, nchan)
                my_freq_array = my_spectra.get_freqarray()
                max_spectrum = np.nanmax(my_dynspec_data.squeeze(), axis=1)
                mean_spectrum = np.nanmean(my_dynspec_data.squeeze(), axis=1)
                median_spectrum = np.nanmedian(my_dynspec_data.squeeze(), axis=1)

                short_file_name = file[:15] + "..." + file[-15:]
                axs.plot(my_freq_array, median_spectrum, label=short_file_name + ' median')
                axs.plot(my_freq_array, mean_spectrum, label=short_file_name + ' mean')
                axs.plot(my_freq_array, max_spectrum, label=short_file_name + ' max')
                axs.legend(loc='upper right')
                axs.set_title(short_file_name)
            else:
                my_spectra.rm_pfb()
                my_dynspec_data = my_spectra.get_data() # (nsubint, nchan)
                my_freq_array = my_spectra.get_freqarray()
                max_spectrum = np.nanmax(my_dynspec_data.squeeze(), axis=1)
                mean_spectrum = np.nanmean(my_dynspec_data.squeeze(), axis=1)
                median_spectrum = np.nanmedian(my_dynspec_data.squeeze(), axis=1)

                short_file_name = file[:15] + "..." + file[-15:]
                axs[index].plot(my_freq_array, median_spectrum, label=short_file_name + ' median')
                axs[index].plot(my_freq_array, mean_spectrum, label=short_file_name + ' mean')
                axs[index].plot(my_freq_array, max_spectrum, label=short_file_name + ' max')
                axs[index].legend(loc='upper right')
                axs[index].set_title(short_file_name)

            

        date = my_spectra.get_start()
        # print(date)
        
        my_dynspec_data = my_spectra.get_data()
        my_freq_array = my_spectra.get_freqarray()
        timearray = my_spectra.get_timearray()
        # print(timearray.iso)
    if len(args.files) == 0:
        pass
    else:
        plt.xlabel("Frequency (MHz)")
        plt.ylabel("Intensity")
        plt.show()
