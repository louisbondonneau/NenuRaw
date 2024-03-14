#!/usr/bin/env python3

from NenuRaw_module import Dynspec

if __name__ == "__main__":
    # wavfile (GUPPI or RAWTF format)
    files = ['/databf/nenufar-pulsar/ES03/2020/05/B0834+06_D20200510T1700_58979_250507_0069_BEAM1.0000.raw']

    # initialisation of the Dynspec object containing the methodes
    my_spectra = Dynspec(files)

    my_spectra.print_info()
