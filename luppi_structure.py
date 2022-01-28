
import numpy as np


class LUPPI:
    def __init__(self, fn_raw, verbose=False):
        """Parse raw file
        """
        self.verbose = verbose

        header_test = [self.__try_header0, self.__try_header1, self.__try_header2, self.__try_header3, self.__try_header4]
        for header_func in header_test:
            result, self.dt_header = header_func(fn_raw)
            if (result):
                self.header = result
                break

        if not (result):
            print('ERROR: the corresponding header is not in raw_structure.py')
            self.print_header(fn_raw)
            exit(0)

        nchan = int(self.header['OBSNCHAN'][9:])
        npol = int(self.header['NPOL'][9:])
        nbit = int(self.header['NBITS'][9:])
        blocksize = int(self.header['BLOCSIZE'][9:])
        overlap = int(self.header['OVERLAP'][9:])  # FFTLEN
        bytespersample = (nbit / 8) * npol * nchan

        if (nbit == 8):
            bit_mode = 'int8'
        elif (nbit == 16):
            bit_mode = 'int16'
        else:
            print("Unknown nbit:")

        print(bit_mode, blocksize)

        # dt_block = np.dtype([('data', bit_mode, int((blocksize-overlap)/(nbit/8))),
        #                     ('overlap', bit_mode, int((overlap)/(nbit/8))),
        #                       ])
        self.dt_block = np.dtype([('header', self.dt_header),
                                  ('data', bit_mode, int(blocksize)),
                                  #('data', bit_mode, int(blocksize-overlap)),
                                  #('overlap', bit_mode, int(overlap)),
                                  ])

        if self.verbose:
            print("  blocksize_header  is: %s -> %.3f sec" % (blocksize, float((blocksize / nchan / npol) - overlap) * 5.12e-6))
            print("  block_size_file   is: %s" % self.dt_block.itemsize)
            print("  overlap    is: %s" % overlap)
            print("  nbit       is: %s" % nbit)
            print("  head_size  is: %s" % self.dt_header.itemsize)

        self.data = self.__open_raw(fn_raw, self.dt_header, self.dt_block)

        self.nchan_file = int(self.header['OBSNCHAN'][9:])
        self.dm = float(self.header['DM'][9:])
        self.nof_blocks = int(self.data['data'].shape[0])
        self.overlap = int(self.header['OVERLAP'][9:])
        self.nof_polcpx = int(self.header['NPOL'][9:])

        self.nschan_file = int(self.data['data'].shape[1] / self.nchan_file)
        self.nschan_file = int(self.nschan_file - self.overlap * self.nof_polcpx)

        self.nschan = self.nschan_file
        self.tbin = float(self.header['TBIN'][9:])
        self.time_int = float(int(self.header['BLOCSIZE'][9:]) / self.nchan_file / self.nof_polcpx) * self.tbin
        self.beamletmin_file = int(self.header['LOWCHAN'][9:])
        self.beamletmax_file = int(self.header['HIGHCHAN'][9:])
        self.chan_bw = float(self.header['CHAN_BW'][9:])
        self.imjd = float(self.header['STT_IMJD'][9:])
        self.smjd = float(self.header['STT_SMJD'][9:])
        self.offmjd = float(self.header['STT_OFFS'][9:])

        # self.patidx = int(self.data[file_num][block_num]['header']['PKTIDX'][9:]) #1585112455.0001945
    def get_patidx(self, data, file_num, block_num):
        return 8 * int(data[file_num][block_num]['header']['PKTIDX'][9:])

    def get_block(self, data, file_num, block_num):
        return data[file_num]['data'][block_num]  # shape is 1D (nchan, nschan, npol)

    def __open_raw(self, fn_raw, dt_header, dt_block):
        """Return a memory mapped file array"""
        with fn_raw.open('rb') as fd_raw:
            fd_raw.seek(0, 2)
            file_size = fd_raw.tell()
            nblock = int(np.floor(file_size / dt_block.itemsize))
            print("  nblock     is: %s" % nblock)
            if((fd_raw.tell() % float(dt_block.itemsize)) != 0):
                print('WARNING: last block is corrupted %f' % (fd_raw.tell() / float(dt_block.itemsize)))
            # print(float(fd_raw.tell())/1024/1024/1024)
        # print(np.__file__)
        data = np.memmap(fn_raw,
                         dtype=dt_block,
                         mode='r',
                         shape=(nblock,)
                         )
        # Reinterpret an integer number of dt_block-length samples as dt_block
        # nof_blocks = tmp.size * tmp.itemsize // (dt_block.itemsize)
        # data = tmp[: nof_blocks * dt_block.itemsize].view(dt_block)
        # print(np.shape(tmp))
        # print(tmp[118]['header'])
        # print(np.shape(tmp['data']))
        return data

    def __try_header0(self, fn_raw):
        dt_header_tmp = np.dtype([('RUN', 'S80'),
                                  ('TELESCOP', 'S80'),
                                  ('FRONTEND', 'S80'),
                                  ('BACKEND', 'S80'),
                                  ('PKTFMT', 'S80'),
                                  ('PKTSIZE', 'S80'),
                                  ('DATAHOST', 'S80'),
                                  ('DATAPORT', 'S80'),
                                  ('NRCVR', 'S80'),
                                  ('NBITS', 'S80'),
                                  ('NPOL', 'S80'),
                                  ('POL_TYPE', 'S80'),
                                  ('FD_POLN', 'S80'),
                                  ('TRK_MODE', 'S80'),
                                  ('N_DS', 'S80'),
                                  ('N_GPU', 'S80'),
                                  ('SRC_NAME', 'S80'),
                                  ('BEAM_NUM', 'S80'),
                                  ('MJDSTART', 'S80'),
                                  ('RA_STR', 'S80'),
                                  ('DEC_STR', 'S80'),
                                  ('DM', 'S80'),
                                  ('RM', 'S80'),
                                  ('OVERLAP', 'S80'),
                                  ('PARFILE', 'S80'),
                                  ('ONLY_I', 'S80'),
                                  ('DS_TIME', 'S80'),
                                  ('DS_FREQ', 'S80'),
                                  ('CHAN_DM', 'S80'),
                                  ('OBS_MODE', 'S80'),
                                  ('TFOLD', 'S80'),
                                  ('OBS_LEN', 'S80'),
                                  ('SCANLEN', 'S80'),
                                  ('TOTNCHAN', 'S80'),
                                  ('OBSBW', 'S80'),
                                  ('OBSNCHAN', 'S80'),
                                  ('LOWCHAN', 'S80'),
                                  ('HIGHCHAN', 'S80'),
                                  ('OBSERVER', 'S80'),
                                  ('PROJID', 'S80'),
                                  ('CAL_MODE', 'S80'),
                                  ('CAL_FREQ', 'S80'),
                                  ('CAL_DCYC', 'S80'),
                                  ('CAL_PHS', 'S80'),
                                  ('STT_IMJD', 'S80'),
                                  ('STT_JDAY', 'S80'),
                                  ('STT_SMJD', 'S80'),
                                  ('STT_OFFS', 'S80'),
                                  ('OFFSET0', 'S80'),
                                  ('SCALE0', 'S80'),
                                  ('OFFSET1', 'S80'),
                                  ('SCALE1', 'S80'),
                                  ('OFFSET2', 'S80'),
                                  ('SCALE2', 'S80'),
                                  ('OFFSET3', 'S80'),
                                  ('SCALE3', 'S80'),
                                  ('ACC_LEN', 'S80'),
                                  ('NBITSADC', 'S80'),
                                  ('PFB_OVER', 'S80'),
                                  ('AZ', 'S80'),
                                  ('ZA', 'S80'),
                                  ('BMAJ', 'S80'),
                                  ('BMIN', 'S80'),
                                  ('LST', 'S80'),
                                  ('TSUBINT', 'S80'),
                                  ('OFFS_SUB', 'S80'),
                                  ('NPOLYCO', 'S80'),
                                  ('NPKT', 'S80'),
                                  ('EQUINOX', 'S80'),
                                  ('OBSFREQ', 'S80'),
                                  ('CHAN_BW', 'S80'),
                                  ('NBIN', 'S80'),
                                  ('SNGPULSE', 'S80'),
                                  ('DATADIR', 'S80'),
                                  ('RAWHOST', 'S80'),
                                  ('DST_IP', 'S80'),
                                  ('RCV_IP', 'S80'),
                                  ('RA', 'S80'),
                                  ('DEC', 'S80'),
                                  ('TBIN', 'S80'),
                                  ('FFTLEN', 'S80'),
                                  ('BLOCSIZE', 'S80'),
                                  ('SCANNUM', 'S80'),
                                  ('BASENAME', 'S80'),
                                  ('DISKSTAT', 'S80'),
                                  ('NETSTAT', 'S80'),
                                  ('DROPAVG', 'S80'),
                                  ('DROPTOT', 'S80'),
                                  ('DROPBLK', 'S80'),
                                  ('STTVALID', 'S80'),
                                  ('DISPSTAT', 'S80'),
                                  ('CURBLOCK', 'S80'),
                                  ('PKTIDX', 'S80'),
                                  ('NDROP', 'S80'),
                                  ('DATATYPE', 'S80'),
                                  ('END', 'S80'),
                                  ])
        print(fn_raw)
        with fn_raw.open('rb') as fd_raw:
            header_tmp = np.frombuffer(fd_raw.read(dt_header_tmp.itemsize),
                                       count=1,
                                       dtype=dt_header_tmp,
                                       )[0]
        c1 = (header_tmp['BLOCSIZE'][:8] == b'BLOCSIZE')
        c2 = (header_tmp['OVERLAP'][:7] == b'OVERLAP')
        c3 = (header_tmp['NBITS'][:5] == b'NBITS')
        c4 = (header_tmp['NPOL'][:4] == b'NPOL')
        c5 = (header_tmp['OBSNCHAN'][:8] == b'OBSNCHAN')
        c6 = (header_tmp['DM'][:2] == b'DM')
        c7 = (header_tmp['LOWCHAN'][:7] == b'LOWCHAN')
        c8 = (header_tmp['HIGHCHAN'][:8] == b'HIGHCHAN')
        c9 = (header_tmp['PKTIDX'][:6] == b'PKTIDX')
        # print(c1, c2, c3, c4, c5, c6, c7, c8)
        if c1 and c2 and c3 and c4 and c5 and c6 and c7 and c8 and c9:
            return header_tmp, dt_header_tmp
        else:
            return False, dt_header_tmp

    def __try_header1(self, fn_raw):
        dt_header_tmp = np.dtype([('RUN', 'S80'),
                                  ('TELESCOP', 'S80'),
                                  ('FRONTEND', 'S80'),
                                  ('BACKEND', 'S80'),
                                  ('PKTFMT', 'S80'),
                                  ('PKTSIZE', 'S80'),
                                  ('DATAHOST', 'S80'),
                                  ('DATAPORT', 'S80'),
                                  ('NRCVR', 'S80'),
                                  ('NBITS', 'S80'),
                                  ('NPOL', 'S80'),
                                  ('POL_TYPE', 'S80'),
                                  ('FD_POLN', 'S80'),
                                  ('TRK_MODE', 'S80'),
                                  ('N_DS', 'S80'),
                                  ('N_GPU', 'S80'),
                                  ('SRC_NAME', 'S80'),
                                  ('BEAM_NUM', 'S80'),
                                  ('MJDSTART', 'S80'),
                                  ('RA_STR', 'S80'),
                                  ('DEC_STR', 'S80'),
                                  ('DM', 'S80'),
                                  ('RM', 'S80'),
                                  ('OVERLAP', 'S80'),
                                  ('PARFILE', 'S80'),
                                  ('ONLY_I', 'S80'),
                                  ('DS_TIME', 'S80'),
                                  ('DS_FREQ', 'S80'),
                                  ('CHAN_DM', 'S80'),
                                  ('OBS_MODE', 'S80'),
                                  ('TFOLD', 'S80'),
                                  ('OBS_LEN', 'S80'),
                                  ('SCANLEN', 'S80'),
                                  ('TOTNCHAN', 'S80'),
                                  ('OBSBW', 'S80'),
                                  ('OBSNCHAN', 'S80'),
                                  ('LOWCHAN', 'S80'),
                                  ('HIGHCHAN', 'S80'),
                                  ('OBSERVER', 'S80'),
                                  ('PROJID', 'S80'),
                                  ('CAL_MODE', 'S80'),
                                  ('CAL_FREQ', 'S80'),
                                  ('CAL_DCYC', 'S80'),
                                  ('CAL_PHS', 'S80'),
                                  ('STT_IMJD', 'S80'),
                                  ('STT_JDAY', 'S80'),
                                  ('STT_SMJD', 'S80'),
                                  ('STT_OFFS', 'S80'),
                                  ('OFFSET0', 'S80'),
                                  ('SCALE0', 'S80'),
                                  ('OFFSET1', 'S80'),
                                  ('SCALE1', 'S80'),
                                  ('OFFSET2', 'S80'),
                                  ('SCALE2', 'S80'),
                                  ('OFFSET3', 'S80'),
                                  ('SCALE3', 'S80'),
                                  ('ACC_LEN', 'S80'),
                                  ('NBITSADC', 'S80'),
                                  ('PFB_OVER', 'S80'),
                                  ('AZ', 'S80'),
                                  ('ZA', 'S80'),
                                  ('BMAJ', 'S80'),
                                  ('BMIN', 'S80'),
                                  ('LST', 'S80'),
                                  ('TSUBINT', 'S80'),
                                  ('OFFS_SUB', 'S80'),
                                  ('NPOLYCO', 'S80'),
                                  ('NPKT', 'S80'),
                                  ('EQUINOX', 'S80'),
                                  ('OBSFREQ', 'S80'),
                                  ('CHAN_BW', 'S80'),
                                  ('NBIN', 'S80'),
                                  ('SNGPULSE', 'S80'),
                                  ('DATADIR', 'S80'),
                                  ('RAWHOST', 'S80'),
                                  ('RAWPORT', 'S80'),
                                  ('RA', 'S80'),
                                  ('DEC', 'S80'),
                                  ('TBIN', 'S80'),
                                  ('FFTLEN', 'S80'),
                                  ('BLOCSIZE', 'S80'),
                                  ('SCANNUM', 'S80'),
                                  ('BASENAME', 'S80'),
                                  ('DISKSTAT', 'S80'),
                                  ('NETSTAT', 'S80'),
                                  ('DROPAVG', 'S80'),
                                  ('DROPTOT', 'S80'),
                                  ('DROPBLK', 'S80'),
                                  ('STTVALID', 'S80'),
                                  ('DISPSTAT', 'S80'),
                                  ('CURBLOCK', 'S80'),
                                  ('PKTIDX', 'S80'),
                                  ('NDROP', 'S80'),
                                  ('DATATYPE', 'S80'),
                                  ('END', 'S80'),
                                  #('data', 'uint8', int(100663296)),
                                  #('overlap', 'uint8', int(100663296)),
                                  #('header2', 'S7600')
                                  ])

        with fn_raw.open('rb') as fd_raw:
            header_tmp = np.frombuffer(fd_raw.read(dt_header_tmp.itemsize),
                                       count=1,
                                       dtype=dt_header_tmp,
                                       )[0]
        c1 = (header_tmp['BLOCSIZE'][:8] == b'BLOCSIZE')
        c2 = (header_tmp['OVERLAP'][:7] == b'OVERLAP')
        c3 = (header_tmp['NBITS'][:5] == b'NBITS')
        c4 = (header_tmp['NPOL'][:4] == b'NPOL')
        c5 = (header_tmp['OBSNCHAN'][:8] == b'OBSNCHAN')
        c6 = (header_tmp['DM'][:2] == b'DM')
        c7 = (header_tmp['LOWCHAN'][:7] == b'LOWCHAN')
        c8 = (header_tmp['HIGHCHAN'][:8] == b'HIGHCHAN')
        c9 = (header_tmp['PKTIDX'][:6] == b'PKTIDX')
        #print(c1, c2, c3, c4, c5, c6, c7, c8)
        if c1 and c2 and c3 and c4 and c5 and c6 and c7 and c8 and c9:
            return header_tmp, dt_header_tmp
        else:
            return False, dt_header_tmp

    def __try_header2(self, fn_raw):
        dt_header_tmp = np.dtype([('RUN', 'S80'),
                                  ('TELESCOP', 'S80'),
                                  ('FRONTEND', 'S80'),
                                  ('BACKEND', 'S80'),
                                  ('PKTFMT', 'S80'),
                                  ('PKTSIZE', 'S80'),
                                  ('DATAHOST', 'S80'),
                                  ('DATAPORT', 'S80'),
                                  ('NRCVR', 'S80'),
                                  ('NBITS', 'S80'),
                                  ('NPOL', 'S80'),
                                  ('POL_TYPE', 'S80'),
                                  ('FD_POLN', 'S80'),
                                  ('TRK_MODE', 'S80'),
                                  ('N_DS', 'S80'),
                                  ('N_GPU', 'S80'),
                                  ('SRC_NAME', 'S80'),
                                  ('BEAM_NUM', 'S80'),
                                  ('MJDSTART', 'S80'),
                                  ('RA_STR', 'S80'),
                                  ('DEC_STR', 'S80'),
                                  ('DM', 'S80'),
                                  ('OVERLAP', 'S80'),
                                  ('PARFILE', 'S80'),
                                  ('ONLY_I', 'S80'),
                                  ('DS_TIME', 'S80'),
                                  ('DS_FREQ', 'S80'),
                                  ('CHAN_DM', 'S80'),
                                  ('OBS_MODE', 'S80'),
                                  ('TFOLD', 'S80'),
                                  ('OBS_LEN', 'S80'),
                                  ('SCANLEN', 'S80'),
                                  ('TOTNCHAN', 'S80'),
                                  ('OBSBW', 'S80'),
                                  ('OBSNCHAN', 'S80'),
                                  ('LOWCHAN', 'S80'),
                                  ('HIGHCHAN', 'S80'),
                                  ('OBSERVER', 'S80'),
                                  ('CAL_MODE', 'S80'),
                                  ('CAL_FREQ', 'S80'),
                                  ('CAL_DCYC', 'S80'),
                                  ('CAL_PHS', 'S80'),
                                  ('STT_IMJD', 'S80'),
                                  ('STT_JDAY', 'S80'),
                                  ('STT_SMJD', 'S80'),
                                  ('STT_OFFS', 'S80'),
                                  ('OFFSET0', 'S80'),
                                  ('SCALE0', 'S80'),
                                  ('OFFSET1', 'S80'),
                                  ('SCALE1', 'S80'),
                                  ('OFFSET2', 'S80'),
                                  ('SCALE2', 'S80'),
                                  ('OFFSET3', 'S80'),
                                  ('SCALE3', 'S80'),
                                  ('ACC_LEN', 'S80'),
                                  ('NBITSADC', 'S80'),
                                  ('PFB_OVER', 'S80'),
                                  ('AZ', 'S80'),
                                  ('ZA', 'S80'),
                                  ('BMAJ', 'S80'),
                                  ('BMIN', 'S80'),
                                  ('LST', 'S80'),
                                  ('TSUBINT', 'S80'),
                                  ('OFFS_SUB', 'S80'),
                                  ('NPOLYCO', 'S80'),
                                  ('NPKT', 'S80'),
                                  ('EQUINOX', 'S80'),
                                  ('OBSFREQ', 'S80'),
                                  ('CHAN_BW', 'S80'),
                                  ('NBIN', 'S80'),
                                  ('SNGPULSE', 'S80'),
                                  ('PROJID', 'S80'),
                                  ('DATADIR', 'S80'),
                                  ('RAWHOST', 'S80'),
                                  ('RAWPORT', 'S80'),
                                  ('RA', 'S80'),
                                  ('DEC', 'S80'),
                                  ('TBIN', 'S80'),
                                  ('FFTLEN', 'S80'),
                                  ('BLOCSIZE', 'S80'),
                                  ('SCANNUM', 'S80'),
                                  ('BASENAME', 'S80'),
                                  ('DISKSTAT', 'S80'),
                                  ('NETSTAT', 'S80'),
                                  ('DROPAVG', 'S80'),
                                  ('DROPTOT', 'S80'),
                                  ('DROPBLK', 'S80'),
                                  ('STTVALID', 'S80'),
                                  ('DISPSTAT', 'S80'),
                                  ('CURBLOCK', 'S80'),
                                  ('PKTIDX', 'S80'),
                                  ('RM', 'S80'),
                                  ('NDROP', 'S80'),
                                  ('DATATYPE', 'S80'),
                                  ('END', 'S80'),
                                  ])

        with fn_raw.open('rb') as fd_raw:
            header_tmp = np.frombuffer(fd_raw.read(dt_header_tmp.itemsize),
                                       count=1,
                                       dtype=dt_header_tmp,
                                       )[0]
        c1 = (header_tmp['BLOCSIZE'][:8] == b'BLOCSIZE')
        c2 = (header_tmp['OVERLAP'][:7] == b'OVERLAP')
        c3 = (header_tmp['NBITS'][:5] == b'NBITS')
        c4 = (header_tmp['NPOL'][:4] == b'NPOL')
        c5 = (header_tmp['OBSNCHAN'][:8] == b'OBSNCHAN')
        c6 = (header_tmp['DM'][:2] == b'DM')
        c7 = (header_tmp['LOWCHAN'][:7] == b'LOWCHAN')
        c8 = (header_tmp['HIGHCHAN'][:8] == b'HIGHCHAN')
        c9 = (header_tmp['PKTIDX'][:6] == b'PKTIDX')
        #print(c1, c2, c3, c4, c5, c6, c7, c8)
        if c1 and c2 and c3 and c4 and c5 and c6 and c7 and c8 and c9:
            return header_tmp, dt_header_tmp
        else:
            return False, dt_header_tmp

    def __try_header3(self, fn_raw):
        dt_header_tmp = np.dtype([('TELESCOP', 'S80'),
                                  ('FRONTEND', 'S80'),
                                  ('BACKEND', 'S80'),
                                  ('PKTFMT', 'S80'),
                                  ('PKTSIZE', 'S80'),
                                  ('DATAHOST', 'S80'),
                                  ('DATAPORT', 'S80'),
                                  ('NRCVR', 'S80'),
                                  ('NBITS', 'S80'),
                                  ('NPOL', 'S80'),
                                  ('POL_TYPE', 'S80'),
                                  ('FD_POLN', 'S80'),
                                  ('TRK_MODE', 'S80'),
                                  ('N_DS', 'S80'),
                                  ('N_GPU', 'S80'),
                                  ('SRC_NAME', 'S80'),
                                  ('BEAM_NUM', 'S80'),
                                  ('MJDSTART', 'S80'),
                                  ('RA_STR', 'S80'),
                                  ('DEC_STR', 'S80'),
                                  ('DM', 'S80'),
                                  ('RM', 'S80'),
                                  ('OVERLAP', 'S80'),
                                  ('PARFILE', 'S80'),
                                  ('ONLY_I', 'S80'),
                                  ('DS_TIME', 'S80'),
                                  ('DS_FREQ', 'S80'),
                                  ('CHAN_DM', 'S80'),
                                  ('OBS_MODE', 'S80'),
                                  ('TFOLD', 'S80'),
                                  ('OBS_LEN', 'S80'),
                                  ('SCANLEN', 'S80'),
                                  ('TOTNCHAN', 'S80'),
                                  ('OBSBW', 'S80'),
                                  ('OBSNCHAN', 'S80'),
                                  ('LOWCHAN', 'S80'),
                                  ('HIGHCHAN', 'S80'),
                                  ('OBSERVER', 'S80'),
                                  ('PROJID', 'S80'),
                                  ('CAL_MODE', 'S80'),
                                  ('CAL_FREQ', 'S80'),
                                  ('CAL_DCYC', 'S80'),
                                  ('CAL_PHS', 'S80'),
                                  ('STT_IMJD', 'S80'),
                                  ('STT_JDAY', 'S80'),
                                  ('STT_SMJD', 'S80'),
                                  ('STT_OFFS', 'S80'),
                                  ('OFFSET0', 'S80'),
                                  ('SCALE0', 'S80'),
                                  ('OFFSET1', 'S80'),
                                  ('SCALE1', 'S80'),
                                  ('OFFSET2', 'S80'),
                                  ('SCALE2', 'S80'),
                                  ('OFFSET3', 'S80'),
                                  ('SCALE3', 'S80'),
                                  ('ACC_LEN', 'S80'),
                                  ('NBITSADC', 'S80'),
                                  ('PFB_OVER', 'S80'),
                                  ('AZ', 'S80'),
                                  ('ZA', 'S80'),
                                  ('BMAJ', 'S80'),
                                  ('BMIN', 'S80'),
                                  ('LST', 'S80'),
                                  ('TSUBINT', 'S80'),
                                  ('OFFS_SUB', 'S80'),
                                  ('NPOLYCO', 'S80'),
                                  ('NPKT', 'S80'),
                                  ('EQUINOX', 'S80'),
                                  ('OBSFREQ', 'S80'),
                                  ('CHAN_BW', 'S80'),
                                  ('NBIN', 'S80'),
                                  ('SNGPULSE', 'S80'),
                                  ('DATADIR', 'S80'),
                                  ('RAWHOST', 'S80'),
                                  ('DST_IP', 'S80'),
                                  ('RCV_IP', 'S80'),
                                  ('RA', 'S80'),
                                  ('DEC', 'S80'),
                                  ('TBIN', 'S80'),
                                  ('FFTLEN', 'S80'),
                                  ('BLOCSIZE', 'S80'),
                                  ('SCANNUM', 'S80'),
                                  ('BASENAME', 'S80'),
                                  ('RUN', 'S80'),
                                  ('DISKSTAT', 'S80'),
                                  ('NETSTAT', 'S80'),
                                  ('DROPAVG', 'S80'),
                                  ('DROPTOT', 'S80'),
                                  ('DROPBLK', 'S80'),
                                  ('STTVALID', 'S80'),
                                  ('DISPSTAT', 'S80'),
                                  ('CURBLOCK', 'S80'),
                                  ('NULLSTAT', 'S80'),
                                  ('PKTIDX', 'S80'),
                                  ('NDROP', 'S80'),
                                  ('DATATYPE', 'S82'),
                                  ('END', 'S80'),
                                  ])

        with fn_raw.open('rb') as fd_raw:
            header_tmp = np.frombuffer(fd_raw.read(dt_header_tmp.itemsize),
                                       count=1,
                                       dtype=dt_header_tmp,
                                       )[0]
        c1 = (header_tmp['BLOCSIZE'][:8] == b'BLOCSIZE')
        c2 = (header_tmp['OVERLAP'][:7] == b'OVERLAP')
        c3 = (header_tmp['NBITS'][:5] == b'NBITS')
        c4 = (header_tmp['NPOL'][:4] == b'NPOL')
        c5 = (header_tmp['OBSNCHAN'][:8] == b'OBSNCHAN')
        c6 = (header_tmp['DM'][:2] == b'DM')
        c7 = (header_tmp['LOWCHAN'][:7] == b'LOWCHAN')
        c8 = (header_tmp['HIGHCHAN'][:8] == b'HIGHCHAN')
        c9 = (header_tmp['PKTIDX'][:6] == b'PKTIDX')
        #print(c1, c2, c3, c4, c5, c6, c7, c8)
        if c1 and c2 and c3 and c4 and c5 and c6 and c7 and c8 and c9:
            return header_tmp, dt_header_tmp
        else:
            return False, dt_header_tmp

    def __try_header4(self, fn_raw):
        dt_header_tmp = np.dtype([('TELESCOP', 'S80'),
                                  ('FRONTEND', 'S80'),
                                  ('BACKEND', 'S80'),
                                  ('PKTFMT', 'S80'),
                                  ('PKTSIZE', 'S80'),
                                  ('DATAHOST', 'S80'),
                                  ('DATAPORT', 'S80'),
                                  ('NRCVR', 'S80'),
                                  ('NBITS', 'S80'),
                                  ('NPOL', 'S80'),
                                  ('POL_TYPE', 'S80'),
                                  ('FD_POLN', 'S80'),
                                  ('TRK_MODE', 'S80'),
                                  ('N_DS', 'S80'),
                                  ('N_GPU', 'S80'),
                                  ('SRC_NAME', 'S80'),
                                  ('BEAM_NUM', 'S80'),
                                  ('MJDSTART', 'S80'),
                                  ('RA_STR', 'S80'),
                                  ('DEC_STR', 'S80'),
                                  ('DM', 'S80'),
                                  ('RM', 'S80'),
                                  ('OVERLAP', 'S80'),
                                  ('PARFILE', 'S80'),
                                  ('ONLY_I', 'S80'),
                                  ('DS_TIME', 'S80'),
                                  ('DS_FREQ', 'S80'),
                                  ('CHAN_DM', 'S80'),
                                  ('OBS_MODE', 'S80'),
                                  ('TFOLD', 'S80'),
                                  ('OBS_LEN', 'S80'),
                                  ('SCANLEN', 'S80'),
                                  ('TOTNCHAN', 'S80'),
                                  ('OBSBW', 'S80'),
                                  ('OBSNCHAN', 'S80'),
                                  ('LOWCHAN', 'S80'),
                                  ('HIGHCHAN', 'S80'),
                                  ('OBSERVER', 'S80'),
                                  ('PROJID', 'S80'),
                                  ('CAL_MODE', 'S80'),
                                  ('CAL_FREQ', 'S80'),
                                  ('CAL_DCYC', 'S80'),
                                  ('CAL_PHS', 'S80'),
                                  ('STT_IMJD', 'S80'),
                                  ('STT_JDAY', 'S80'),
                                  ('STT_SMJD', 'S80'),
                                  ('STT_OFFS', 'S80'),
                                  ('OFFSET0', 'S80'),
                                  ('SCALE0', 'S80'),
                                  ('OFFSET1', 'S80'),
                                  ('SCALE1', 'S80'),
                                  ('OFFSET2', 'S80'),
                                  ('SCALE2', 'S80'),
                                  ('OFFSET3', 'S80'),
                                  ('SCALE3', 'S80'),
                                  ('ACC_LEN', 'S80'),
                                  ('NBITSADC', 'S80'),
                                  ('PFB_OVER', 'S80'),
                                  ('AZ', 'S80'),
                                  ('ZA', 'S80'),
                                  ('BMAJ', 'S80'),
                                  ('BMIN', 'S80'),
                                  ('LST', 'S80'),
                                  ('TSUBINT', 'S80'),
                                  ('OFFS_SUB', 'S80'),
                                  ('NPOLYCO', 'S80'),
                                  ('NPKT', 'S80'),
                                  ('EQUINOX', 'S80'),
                                  ('OBSFREQ', 'S80'),
                                  ('CHAN_BW', 'S80'),
                                  ('NBIN', 'S80'),
                                  ('SNGPULSE', 'S80'),
                                  ('DATADIR', 'S80'),
                                  ('RAWHOST', 'S80'),
                                  ('DST_IP', 'S80'),
                                  ('RCV_IP', 'S80'),
                                  ('RA', 'S80'),
                                  ('DEC', 'S80'),
                                  ('TBIN', 'S80'),
                                  ('FFTLEN', 'S80'),
                                  ('BLOCSIZE', 'S80'),
                                  ('SCANNUM', 'S80'),
                                  ('BASENAME', 'S80'),
                                  ('RUN', 'S80'),
                                  ('DISKSTAT', 'S80'),
                                  ('NETSTAT', 'S80'),
                                  ('DROPAVG', 'S80'),
                                  ('DROPTOT', 'S80'),
                                  ('DROPBLK', 'S80'),
                                  ('STTVALID', 'S80'),
                                  ('DISPSTAT', 'S80'),
                                  ('CURBLOCK', 'S80'),
                                  ('PKTIDX', 'S80'),
                                  ('NDROP', 'S80'),
                                  ('DATATYPE', 'S80'),
                                  ('END', 'S80'),
                                  ])

        with fn_raw.open('rb') as fd_raw:
            header_tmp = np.frombuffer(fd_raw.read(dt_header_tmp.itemsize),
                                       count=1,
                                       dtype=dt_header_tmp,
                                       )[0]

        c1 = (header_tmp['BLOCSIZE'][:8] == b'BLOCSIZE')
        c2 = (header_tmp['OVERLAP'][:7] == b'OVERLAP')
        c3 = (header_tmp['NBITS'][:5] == b'NBITS')
        c4 = (header_tmp['NPOL'][:4] == b'NPOL')
        c5 = (header_tmp['OBSNCHAN'][:8] == b'OBSNCHAN')
        c6 = (header_tmp['DM'][:2] == b'DM')
        c7 = (header_tmp['LOWCHAN'][:7] == b'LOWCHAN')
        c8 = (header_tmp['HIGHCHAN'][:8] == b'HIGHCHAN')
        c9 = (header_tmp['PKTIDX'][:6] == b'PKTIDX')
        #print(c1, c2, c3, c4, c5, c6, c7, c8)
        if c1 and c2 and c3 and c4 and c5 and c6 and c7 and c8 and c9:
            return header_tmp, dt_header_tmp
        else:
            return False, dt_header_tmp

    def print_header(self, fn_raw):
        dt_header_tmp = np.dtype([('HEADER', 'S9000')])
        with fn_raw.open('rb') as fd_raw:
            header_tmp = np.frombuffer(fd_raw.read(dt_header_tmp.itemsize),
                                       count=1,
                                       dtype=dt_header_tmp,
                                       )[0]
            header_tmp = str(header_tmp).replace('\\', '')
        header_liste = [header_tmp[i:i + 80] for i in range(0, len(header_tmp), 80)]
        # for i in range(int(7600/80)):
        #    print(header_liste[i])
        for i in range(int(9000 / 80)):
            print(header_liste[i])

    def try_luppi(self, fn_raw):
        header_test = [self.__try_header0, self.__try_header1, self.__try_header2, self.__try_header3, self.__try_header4]
        for header_func in header_test:
            result, dt_header = header_func(self, fn_raw)
            if (result):
                return True
        return False
