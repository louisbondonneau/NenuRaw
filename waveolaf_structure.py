
import numpy as np

from datetime import datetime
from astropy.time import Time
import matplotlib.pyplot as plt


class WaveOlaf:
    def __init__(self, fn_raw, verbose=False):
        """Parse raw file
        """
        self.verbose = verbose
        self.header, self.dt_header = self.__header(fn_raw)

        npol = 4
        overlap = 0
        nchan = int(self.header['nof_beamlets'])
        nbit = int(8)
        nschan_file = int(self.header['nofblock']) * npol
        blocksize = nchan * nschan_file

        if (nbit == 8):
            bit_mode = 'int8'
        elif (nbit == 16):
            bit_mode = 'int16'
        else:
            print("Unknown nbit:")

        self.dt_block = np.dtype([('version_id', 'uint8'),
                                  ('source_info', 'uint16'),
                                  ('configuration_id', 'uint8'),
                                  ('station_id', 'uint16'),
                                  ('nof_beamlets', 'uint8'),
                                  ('nofblock', 'uint8'),
                                  ('Timestamp', 'int32'),
                                  ('BSN', 'int32'),
                                  ('data', 'uint8', (int(blocksize * nbit / 8))),
                                  ])

        if self.verbose:
            print("  blocksize_header  is: %s -> %.8f sec" % (blocksize, float((blocksize / float(nchan) / npol) - overlap) * 5.12e-6))
            print("  blocksize_file    is: %s" % self.dt_block['data'].itemsize)
            print("  overlap    is: %s" % overlap)
            print("  nbit       is: %s" % nbit)
            print("  head_size  is: %s" % self.dt_header.itemsize)

        self.data = self.__open_raw(fn_raw, self.dt_header, self.dt_block)
        self.nchan_file = int(self.header['nof_beamlets'])
        self.dm = 0
        self.nof_blocks = int(self.data['data'].shape[0])
        self.nschan = int(self.data['data'].shape[1] / self.nchan_file)
        self.nschan_file = int(self.data['data'].shape[1] / self.nchan_file)
        self.nof_polcpx = int(4)
        self.overlap = 0
        self.tbin = float(5.12e-6)
        self.time_int = float((self.nschan_file / self.nof_polcpx) - self.overlap) * self.tbin
        self.beamletmin_file = int(200)
        self.beamletmax_file = self.beamletmin_file + self.nchan_file - 1
        self.chan_bw = float(200. / 1024)

        date_time = datetime.fromtimestamp(self.header['Timestamp']).isoformat()
        date_time = Time(date_time, format='isot', scale='utc')

        self.imjd = np.floor(date_time.mjd)
        self.smjd = np.round((date_time.mjd % 1) * 3600 * 24)
        if(self.header['Timestamp'] % 1 != 0):
            self.offmjd = self.header['BSN'] * self.tbin + self.tbin / 2
        else:
            self.offmjd = self.header['BSN'] * self.tbin

    def get_patidx(self, data, file_num, block_num):
        ts_0 = data[file_num]['Timestamp'][0]
        if(ts_0 % 1 != 0):
            ts_0 += self.tbin / 2
        ts_1 = data[file_num]['Timestamp'][block_num]
        if(ts_1 % 1 != 0):
            ts_1 += self.tbin / 2

        ts_from_start = ts_1 - ts_0
        return int(ts_from_start / self.tbin + data[file_num]['BSN'][block_num])

    def get_block(self, data, file_num, block_num):
        out = data[file_num]['data'][block_num]
        out.shape = (int(self.nschan_file / self.nof_polcpx), self.nchan_file, self.nof_polcpx)
        out = out.transpose((1, 0, 2))
        return out  # np.ascontiguousarray(out) #shape is 1D (nchan, nschan, npol)

    def __open_raw(self, fn_raw, dt_header, dt_block):
        data = np.memmap(fn_raw.as_posix(),
                         dtype=dt_block,
                         mode='r',
                         offset=0,
                         )
        return data

    def __header(self, fn_raw):
        dt_header_tmp = np.dtype([('version_id', 'uint8'),
                                  ('source_info', 'uint16'),
                                  ('configuration_id', 'uint8'),
                                  ('station_id', 'uint16'),
                                  ('nof_beamlets', 'uint8'),
                                  ('nofblock', 'uint8'),
                                  ('Timestamp', 'int32'),
                                  ('BSN', 'int32'),
                                  ])
        with fn_raw.open('rb') as fd_raw:
            header_tmp = np.frombuffer(fd_raw.read(dt_header_tmp.itemsize),
                                       count=1,
                                       dtype=dt_header_tmp,
                                       )[0]

        return header_tmp, dt_header_tmp

    def try_waveolaf(fn_raw):

        dt_header_tmp = np.dtype([('version_id', 'uint8'),
                                  ('source_info', 'uint16'),
                                  ('configuration_id', 'uint8'),
                                  ('station_id', 'uint16'),
                                  ('nof_beamlets', 'uint8'),
                                  ('nofblock', 'uint8'),
                                  ('Timestamp', 'int32'),
                                  ('BSN', 'int32')
                                  ])
        try:
            with fn_raw.open('rb') as fd_raw:
                header_tmp = np.frombuffer(fd_raw.read(dt_header_tmp.itemsize),
                                           count=1,
                                           dtype=dt_header_tmp,
                                           )[0]
            print(header_tmp['nofblock'])
            print(header_tmp['nof_beamlets'])
            blocksize = 4 * header_tmp['nofblock'] * header_tmp['nof_beamlets']
            dt_header_tmp = np.dtype([('version_id', 'uint8'),
                                      ('source_info', 'uint16'),
                                      ('configuration_id', 'uint8'),
                                      ('station_id', 'uint16'),
                                      ('nof_beamlets', 'uint8'),
                                      ('nofblock', 'uint8'),
                                      ('Timestamp', 'int32'),
                                      ('BSN', 'int32'),
                                      ('data', 'uint8', (int(blocksize))),
                                      ])
            data = np.memmap(fn_raw.as_posix(),
                             dtype=dt_header_tmp,
                             mode='r',
                             offset=0,
                             )
        except ValueError:
            print("waveolaf is not a valid format.  Try again...")
            return False
        return True
