
import numpy as np

import datetime
from astropy.time import Time


class Rawtf:
    def __init__(self, fn_raw, verbose=False):
        """Parse raw file
        """
        self.verbose = verbose

        self.header, self.dt_header = self.__header(fn_raw)

        npol = 4
        overlap = 0
        nchan = int(self.header['OBSNCHAN'])
        nbit = int(8 * self.header['BYTES_SAMPLE'] / npol)
        nschan_file = int(self.header['NSCHAN']) * npol
        blocksize = nchan * nschan_file

        if (nbit == 8):
            bit_mode = 'int8'
        elif (nbit == 16):
            bit_mode = 'int16'
        else:
            print("Unknown nbit:")

        print(bit_mode, blocksize)

        # self.dt_block = np.dtype([('header', self.dt_header),
        #                    ('data', bit_mode, int(blocksize)),
        #                    #('data', bit_mode, int(blocksize-overlap)),
        #                    #('overlap', bit_mode, int(overlap)),
        #                       ])

        self.dt_block = np.dtype([('eisb', 'uint64'),
                                  ('tsb', 'uint64'),
                                  ('bsnb', 'uint64'),
                                  ('data', 'int8', (int(blocksize * nbit / 8),)),
                                  ])

        if self.verbose:
            print("  blocksize_header  is: %s -> %.3f sec" % (blocksize, float((blocksize / nchan / npol) - overlap) * 5.12e-6))
            print("  blocksize_file    is: %s" % self.dt_block['data'].itemsize)
            print("  overlap    is: %s" % overlap)
            print("  nbit       is: %s" % nbit)
            print("  head_size  is: %s" % self.dt_header.itemsize)

        self.data = self.__open_raw(fn_raw, self.dt_header, self.dt_block)
        self.nchan_file = int(self.header['OBSNCHAN'])
        self.dm = 0
        self.nof_blocks = int(self.data['data'].shape[0])
        self.nschan = int(self.data['data'].shape[1] / self.nchan_file)
        self.nschan_file = int(self.data['data'].shape[1] / self.nchan_file)
        self.nof_polcpx = int(self.nschan_file / int(self.header['NSCHAN']))
        self.overlap = 0
        self.tbin = float(5.12e-6)
        self.time_int = float((self.nschan_file / self.nof_polcpx) - self.overlap) * self.tbin
        self.beamletmin_file = int(self.header['lbc_alloc'][0]['chan'])
        self.beamletmax_file = int(self.header['lbc_alloc'][int(self.nchan_file - 1)]['chan'])
        self.chan_bw = float(200. / 1024)

        date_time = datetime.fromtimestamp(self.data['tsb'][0]).isoformat()
        date_time = Time(date_time, format='isot', scale='utc')

        self.imjd = np.floor(date_time.mjd)
        self.smjd = np.round((date_time.mjd % 1) * 3600 * 24)
        self.offmjd = self.data['bsnb'][0] * self.tbin

        # self.patidx = int(self.data[file_num][block_num]['header']['PKTIDX'][9:]) #1585112455.0001945
    def get_patidx(self, data, file_num, block_num):
        return int(data[file_num]['eisb'][block_num])

    def get_block(self, data, file_num, block_num):
        out = data[file_num]['data'][block_num]
        out.shape = (int(self.nschan_file / self.nof_polcpx), self.nchan_file, self.nof_polcpx)
        out = out.transpose((1, 0, 2))
        return np.ascontiguousarray(out)  # shape is 1D (nchan, nschan, npol)

    def __open_raw(self, fn_raw, dt_header, dt_block):
        data = np.memmap(fn_raw.as_posix(),
                         dtype=dt_block,
                         mode='r',
                         offset=dt_header.itemsize,
                         )
        return data

    def __header(self, fn_raw):
        dt_header_tmp = np.dtype([('OBSNCHAN', 'int32'),          # NUMBER_OF_BEAMLET_PER_BANK = number of channels
                                  ('NSCHAN', 'int32'),     # NUMBER of SAMPLES (fftlen*nfft)
                                  ('BYTES_SAMPLE', 'int32'),  # BYTES per SAMPLE (4/8 for 8/16bits data)
                                  ])

        with fn_raw.open('rb') as fd_raw:
            header_tmp = np.frombuffer(fd_raw.read(dt_header_tmp.itemsize),
                                       count=1,
                                       dtype=dt_header_tmp,
                                       )[0]

        nchan = int(header_tmp['OBSNCHAN'])

        dt_lane_beam_chan = np.dtype([('lane', 'int32'),
                                      ('beam', 'int32'),
                                      ('chan', 'int32'),
                                      ])

        dt_header_tmp = np.dtype([('OBSNCHAN', 'int32'),          # NUMBER_OF_BEAMLET_PER_BANK = number of channels
                                  ('NSCHAN', 'int32'),     # NUMBER of SAMPLES (fftlen*nfft)
                                  ('BYTES_SAMPLE', 'int32'),  # BYTES per SAMPLE (4/8 for 8/16bits data)
                                  ('lbc_alloc', dt_lane_beam_chan, (nchan,)),
                                  ])

        with fn_raw.open('rb') as fd_raw:
            header_tmp = np.frombuffer(fd_raw.read(dt_header_tmp.itemsize),
                                       count=1,
                                       dtype=dt_header_tmp,
                                       )[0]

        return header_tmp, dt_header_tmp

    def try_rawtf(fn_raw):
        dt_header_tmp = np.dtype([('OBSNCHAN', 'int32'),          # NUMBER_OF_BEAMLET_PER_BANK = number of channels
                                  ('NSCHAN', 'int32'),     # NUMBER of SAMPLES (fftlen*nfft)
                                  ('BYTES_SAMPLE', 'int32'),  # BYTES per SAMPLE (4/8 for 8/16bits data)
                                  ])
        try:
            with fn_raw.open('rb') as fd_raw:
                header_tmp = np.frombuffer(fd_raw.read(dt_header_tmp.itemsize),
                                           count=1,
                                           dtype=dt_header_tmp,
                                           )[0]
            nchan = int(header_tmp['OBSNCHAN'])
            dt_lane_beam_chan = np.dtype([('lane', 'int32'),
                                          ('beam', 'int32'),
                                          ('chan', 'int32'),
                                          ])
            dt_header_tmp = np.dtype([('OBSNCHAN', 'int32'),          # NUMBER_OF_BEAMLET_PER_BANK = number of channels
                                      ('NSCHAN', 'int32'),     # NUMBER of SAMPLES (fftlen*nfft)
                                      ('BYTES_SAMPLE', 'int32'),  # BYTES per SAMPLE (4/8 for 8/16bits data)
                                      ('lbc_alloc', dt_lane_beam_chan, (nchan,)),
                                      ])
            with fn_raw.open('rb') as fd_raw:
                header_tmp = np.frombuffer(fd_raw.read(dt_header_tmp.itemsize),
                                           count=1,
                                           dtype=dt_header_tmp,
                                           )[0]
        except ValueError:
            print("rawtf is not a valid format.  Try again...")
            return False
        return True
