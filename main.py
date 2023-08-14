#!/usr/bin/env python3

###############################################################################
#
# Copyright (C) 2019
# Station de Radioastronomie de Nançay,
# Observatoire de Paris, PSL Research University, CNRS, Univ. Orléans, OSUC,
# 18330 Nançay, France
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################
#
# Description:
# Lib that can read GUPPI-raw output files
#
###############################################################################

import argparse

from NenuRaw import dynspec_utils
from NenuRaw import wav_utils


parser = argparse.ArgumentParser(description="This code will plot GUPPI raw files")
# parser.add_argument('-gui', dest='gui', action='store_true', default=False,
#                    help="Open the matplotlib graphical user interface")
parser.add_argument('INPUT_ARCHIVE', nargs='+', help="Path to the Archives")

args = parser.parse_args()


if __name__ == "__main__":
    files = args.INPUT_ARCHIVE

    my_spectra = dynspec_utils.Dynspec(files,
                                       verbose=True,
                                       freq_start=99,
                                       freq_end=0,
                                       # start=10,
                                       # end=10.5,
                                       # start=10.1,
                                       # end=10.35
                                       # start=4.6,
                                       # end=4.8,
                                       start=4.7,
                                       end=4.75,
                                       # block_start=0,
                                       # block_end=1,
                                       )

    my_wav_obj = wav_utils.Wav()

    # my_spectra.dm = 12.44
    # my_spectra.new_fourier_methode(my_wav_obj.wav_cleaning_freq)
    my_spectra.new_fourier_methode(my_wav_obj.coherent_dedisp)

    my_spectra.fourier_computation(fftlen=1, ds_ms=0.2, df=1, pol="I")  # I, Q, U, L, V, XX, YY
    # my_spectra.plot_spectra()

    # my_spectra.clean(threshold=5)
    # my_spectra.clean(threshold=10)
    # my_spectra.plot_spectra()
    # my_spectra.plot_dynspec()
    # my_spectra.plot_spectra()
    # my_spectra.clean(threshold=5)
    # my_spectra.plot_spectra()
    # my_spectra.plot_dynspec()

    my_spectra.clean(threshold=10)
    my_spectra.rm_baseline()
    # my_spectra.plot_dynspec()
    my_spectra.clean(threshold=4)
    my_spectra.rm_baseline()
    my_spectra.plot_dynspec()

    # my_spectra.rm_pfb()
    #my_spectra.clean(threshold = 20)
    # for i in range(3):
    #   my_spectra.rm_dm0()
    #   my_spectra.rm_baseline()
    #   my_spectra.clean(threshold = 10)

    # my_spectra.dedisperse()

    # my_spectra.plot_dynspec()
