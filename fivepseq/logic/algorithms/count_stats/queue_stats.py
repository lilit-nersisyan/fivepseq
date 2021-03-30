import logging
from math import floor

import pandas as pd
import numpy as np
from fivepseq.logic.structures.codons import Codons

from fivepseq.logic.algorithms.count_stats.count_stats import CountStats
from fivepseq.logic.structures.fivepseq_counts import FivePSeqCounts, CountManager


class QueueStats:
    fivepseq_counts = None
    fivepseq_out = None
    config = None
    logger = None

    l_coef = 2  # the minimum length of the desired fragments: length = l_coef*period (default: 2)
    wm = 3  # window multiplier for the seeds (default: 3)
    ft = 20  # FFT signal threshold (default: 10)
    tolerance_coef = 0.12  # fraction of allowed variability from desired periodicity (default: 0.12)

    def __init__(self, fivepseq_counts, fivepseq_out, config):
        self.fivepseq_counts = fivepseq_counts
        self.fivepseq_out = fivepseq_out
        self.config = config
        self.logger = logging.getLogger(config.FIVEPSEQ_LOGGER)

    def store_periods(self, period, asite_frag=None):
        """Stores the periods detected by detect_periods() function and stores those in the given dataframe (in-place?).

        :param: period: int a number indicating the length of the periodicity to search for
        :param: asite_frag: int length from A site to ribosome protection site (computed from period if not provided)
        :return: df: pd.DataFrame: a dataframe containing columns ["tnum", "start", "stop", "period", "score", "length", "pos"]
                tnum = transcript number
                start = start of the periodicity window
                stop = stop of the periodicity window
                period = the value of period
                score = the strength of FFT transform signal
                length = length of the window
                pos_s = the position of the window start relative to CDS, in fractions (0 - start, 1 - stop, <0 - 5'UTR, >1 - 3'UTR)
                pos_e = the position of the window end relative to CDS, in fractions (0 - start, 1 - stop, <0 - 5'UTR, >1 - 3'UTR)
                peaks = [start, p1,p2,p3..., stop] the list of the peaks within the queue sequence
                context = period-long sequence from the last peak
                motif = tricodon in the presumed EPA site of the first stalled ribosome
        """

        rq_df = pd.DataFrame(columns=["tnum", "transcript_start", "transcript_stop",
                                      "peak_start", "peak_end",
                                      "peak_genome_start", "peak_genome_end",
                                      "period", "score", "length",
                                      "pos_s", "pos_e",
                                      "from_start", "from_stop", "from_A", "peaks",
                                      "context", "motif", "E", "P", "A"])

        logging.info("Queue analysis started")

        logging.info("Options: ribosome size = %d" %
                     (period))

        count_vector_list = self.fivepseq_counts.get_count_vector_list(FivePSeqCounts.FULL_LENGTH)
        span_size = self.fivepseq_counts.annotation.span_size
        cnt = 0
        ta = self.fivepseq_counts.annotation.get_transcript_assembly(span_size=span_size)
        for ivec in range(len(count_vector_list)):
            vec = count_vector_list[ivec]

            windows = self.get_periods(vec, period)

            if len(windows) > 0:
                # attach periods to df
                for window in windows:
                    vl = len(vec) - span_size * 2
                    s = window[0]
                    e = window[1]
                    peaks = self.get_peaks(vec[s:e], s, e, period)

                    if len(peaks) > 1 or len(peaks) > len(window) / period / 2:

                        real_peaks = [p - span_size for p in peaks]
                        e = peaks[len(peaks) - 1]  # adjust window end to the last peak
                        s = peaks[0]  # adjust window start to the first peak

                        wl = e - s
                        pos_s = float(s - span_size) / vl
                        pos_e = float(e - span_size) / vl
                        from_start = s - span_size
                        from_stop = e - span_size - vl
                        real_s = s - span_size
                        real_e = e - span_size
                        genome_s = ta[ivec].cds_genome_start + real_s
                        genome_e = ta[ivec].cds_genome_start + real_e

                        seq = self.fivepseq_counts.get_sequence(ta[ivec], span_size, span_size)

                        e = peaks[len(peaks) - 1]

                        if e + period <= len(seq):
                            context = seq[e:(e + period)]
                        else:
                            context = seq[e:len(seq)]

                        if asite_frag is None:
                            asite_frag = period - 13  # inferred from 30 and -17 for yeast, but can be changed for other species

                        if e + asite_frag <= len(seq):
                            motif = seq[(e + asite_frag - 6):(e + asite_frag + 3)]  # [+11 : 20) for -17
                        else:
                            motif = ""

                        rq_df = rq_df.append(
                            pd.Series([ivec, ta[ivec].cds_genome_start, ta[ivec].cds_genome_end,
                                       real_s, real_e,
                                       genome_s, genome_e,
                                       window[2], window[3], wl,
                                       pos_s, pos_e, from_start, from_stop, from_stop + 3,
                                       ','.join([str(p) for p in real_peaks]),
                                       context, motif, motif[0:3], motif[3:6], motif[6:9]],
                                      rq_df.columns),
                            ignore_index=True)

            cnt += 1

        logging.info("")
        return (rq_df)

    def get_periods(self, vec, period, quiet=True):
        """ Detects windows in the given vector that contain the desired periodicity.

        :param: vec: [int] a vector of int counts
        :param: period: int a number indicating the length of the periodicity to search for
        :param: l_coef: int indicates the minimum length of the desired fragments: length = l_coef*period (default: 2)
        :param: wm: int window multiplier for the seeds (default: 3)
        :param: ft: int FFT signal threshold (default: 10)
        :param: tolerance_coef: float the fraction of allowed variability from desired periodicity (default: 0.12)

        :return: a list of vectors with fields ["start", "stop", "period", "signal"]

        Constants:
            - wm (window multiplier) - the length of the window is determined by wm*period
            - ft (FFT threshold) - signals greater than this are considered as periodic


        - Divide the whole vector to windows of length 3*period, where 3 is the starting value of wm.
        - On each window:
            -- if the count in the region is greater than or equal to wm-1, check FFT
            -- If FFT is greater than FFT threshold (ft), mark the window for futher processing
        - Collect all the windows with current wm
            -- if the window is isolated, extend its range by period - first from left, then from right - until the strength of the signal continues increasing
            -- if the distance between two adjacent windows is less than or equal to period, merge them
            -- determine the min and max highest points from left and right of each window and slice them
        - Return a list of a list of vectors with fields ["start", "stop", "signal"] for each sliced window


        """
        inc = self.wm * period

        # start and end coordinates of initial seed windows
        s = np.arange(0, len(vec), inc)
        e = s + inc
        e[-1] = len(vec)

        # check periodicity in each window

        valid_windows = []  # vector storing indices of windows with periodicity
        for i in range(len(s)):
            window = vec[s[i]:e[i]]
            if len(window) > period:
                if self.is_valid(window, period, s[i], e[i]):
                    valid_windows.append(i)

            # extend the valid windows
        e_s = []
        e_e = []
        for i in valid_windows:
            window = vec[s[i]:e[i]]
            e_window = self.extend_window(vec, window, s[i], e[i], period, quiet=quiet)
            e_s.append(e_window[0])
            e_e.append(e_window[1])

        # merge extended windows

        m_windows = self.merge_ranges(e_s, e_e, period)
        m_s = m_windows[0]
        m_e = m_windows[1]

        # slice on maximum endpoints
        s_s = []
        s_e = []

        for i in range(len(m_s)):
            s_window = self.slice_range(vec, m_s[i], m_e[i], period)
            s_s.append(s_window[0])
            s_e.append(s_window[1])

        # return slices with at least 2*period length
        result = []
        for i in range(len(s_s)):
            if s_e[i] - s_s[i] >= self.l_coef * period:
                ps = self.get_fft_signal(vec[s_s[i]:s_e[i]], period)
                if ps[0] > self.ft:
                    result.append([s_s[i], s_e[i], ps[1], ps[0]])

        return (result)

    def slice_range(self, vec, s, e, period):
        """Slices the vector so that its endpoints correspond to peak counts.

        :param: vec: [int] a vector of counts
        :param: s: int the start coordinate of the range
        :param: e: int the end coordinate of the range
        :param: period: int the value of periodicity

        :return: [int] the start and end coordinates of the slice
        """

        o_signal = self.get_fft_signal(vec[s:e], period)[0]

        cursor = s

        while True:
            n_signal = self.get_fft_signal(vec[cursor + 1: e], period)[0]
            if n_signal < o_signal:
                break
            o_signal = n_signal
            cursor += 1

        s = cursor
        cursor = e

        while True:
            n_signal = self.get_fft_signal(vec[s: cursor - 1], period)[0]
            if n_signal < o_signal:
                break
            o_signal = n_signal
            cursor -= 1

        e = cursor

        return (s, e)

    def merge_ranges(self, s, e, period):
        """ Given lists of start and end coordinates, returns new lists of starts and ends of ranges, where the distance betweeen two adjacent ranges is at least equal to period.
        :param: s: [int] a list of start coordinates
        :param: e: [int] a list of end coordinates
        :param: period: int minimium distance between ranges

        :return: [[int], [int]] list of novel lists of start and end coordinates for the merged ranges

        """
        if len(s) < 2:
            return (s, e)

        m_s = []
        m_e = []
        new_range = True
        for i in range(len(s) - 1):
            j = i + 1
            if new_range:
                m_s.append(s[i])
            if s[j] - e[i] <= period:
                new_range = False
            else:
                new_range = True
                m_e.append(e[i])
                if j == len(s) - 1:
                    m_s.append(s[j])
        m_e.append(e[j])

        return (m_s, m_e)

    def extend_window(self, vec, window, s, e, period, quiet=True):
        """ Extends the range of current window to increase periodicity signal

        :param: [int] the full vector of counts
        :param: window: [int] vector of counts in the specified window range
        :param: wm: int windows multiplier (explained above)
        :param: period: int the desired periodicity
        :param: ft: thershold for FFT signal
        :param: quiet: logical should FFT stats be not printed

        :return: [int] start and end positions of the new window
        """
        # obtain initial stats
        sp = self.get_fft_signal(window, period)

        o_signal = sp[0]
        o_period = sp[1]

        if not quiet:
            print("Extending window [%d-%d]: with initial signal %f at periodicity %f" % (s, e, o_signal, o_period))

        extend = True
        n_s = s
        o_s = s
        n_e = e
        o_e = e
        while extend:
            o_s = n_s
            n_s = n_s - period
            if n_s <= 0:
                n_s = 0
                extend = False
            n_window = vec[n_s:n_e]
            n_signal = self.get_fft_signal(n_window, period)[0]
            if n_signal < o_signal:
                if not quiet:
                    print("L: Worse signal for [%d-%d]: n(%.2f) versus o(%.2f)" % (n_s, n_e, n_signal, o_signal))
                extend = False
                n_s = o_s
            else:
                if not quiet:
                    print("L: Better signal for [%d-%d]: n(%.2f) versus o(%.2f)" % (n_s, n_e, n_signal, o_signal))
                o_signal = n_signal

        extend = True
        while extend:
            o_e = n_e
            n_e = n_e + period
            if n_e >= len(vec):
                n_e = len(vec)
                extend = False
            n_window = vec[n_s:n_e]
            n_signal = self.get_fft_signal(n_window, period)[0]
            if n_signal < o_signal:
                if not quiet:
                    print("R: Worse signal for [%d-%d]: n(%.2f) versus o(%.2f)" % (n_s, n_e, n_signal, o_signal))
                extend = False
                n_e = o_e
            else:
                if not quiet:
                    print("R: Better signal for [%d-%d]: n(%.2f) versus o(%.2f)" % (n_s, n_e, n_signal, o_signal))
                o_signal = n_signal

        return (n_s, n_e)

    def is_valid(self, window, period, start, stop, quiet=True):
        """ Checks if current window contains counts with given periodicity

        :param: window: [int] vector of counts
        :param: wm: int windows multiplier (explained above)
        :param: period: int the desired periodicity
        :param: ft: threshold for FFT signal
        :param: quiet: logical should FFT stats be not printed

        :return: logical indicating if the window contains desired periodicity
        """
        if len(window) > (self.wm - 1) * period:
            if sum(window) >= self.wm - 1:
                signal = self.get_fft_signal(window, period)[0]

                if signal > self.ft:
                    peaks = self.get_peaks(window, start, stop, period)
                    # number of real peaks should be more than half of possible peaks
                    valid = len(peaks) > len(window) / period / 2
                    return (valid)

        return (False)

    def get_peaks(self, vec, start, stop, period):
        d = np.arange(start, stop)
        count_series = pd.Series(data=vec, index=d)
        peak_pvalues = CountStats.poisson_df_from_count_series(count_series)

        # distances between peaks
        peaks = []
        for d1 in sorted(peak_pvalues.D):
            for d2 in sorted(peak_pvalues.D):
                dif = np.abs(d2 - d1)
                if dif == period or dif == period + 3:
                    if d1 not in peaks:
                        peaks.append(d1)
                    if d2 not in peaks:
                        peaks.append(d2)

        return sorted(peaks)

    def get_fft_signal(self, vec, period, n=3, quiet=True):
        """ Returns the maximum FFT signal observed for the given periodicity within the given vector.

        :param: vec: [int] vector of counts
        :param: period: int the desired periodicity

        :return: float max FFT signal for given periodicity or 0 if periodicity is not observed and corresponding periodicity
        """
        max_signal = 0
        max_ind = -1

        tolerance = self.tolerance_coef * period

        try:

            fft_stats = self.fft_stats_on_vector(vec, 3)
            if not quiet:
                print("Periods: ", fft_stats[1])
                print("Signals: ", fft_stats[2])
                print("Scales:  ", fft_stats[3])

            for j in range(n):
                if abs(fft_stats[1][j] - period) <= tolerance:
                    signal = fft_stats[2][j]
                    if signal > max_signal:
                        max_signal = signal
                        max_ind = j
        except:
            return (0, 0)

        if max_ind == -1:
            max_period = 0
        else:
            max_period = fft_stats[1][max_ind]

        return (max_signal, max_period)

    def fft_stats_on_vector(self, count_vector, n=1):
        # get absolute values of signals for each wave, and take wave-frequencies less than half of the input vector
        fft_abs = np.abs(np.fft.fft(count_vector))
        fft_abs = fft_abs[1:(len(fft_abs) // 2 + 1)]

        # sort the signals from smallest to largest, and take the largest n values
        ind = np.flip(np.argsort(fft_abs)[-1 * n:], 0)
        signals = [fft_abs[i] for i in ind]

        # scale each signal according to the mean of all the signals
        scales = [signal / np.mean(fft_abs) for signal in signals]

        # compute the actual wavelengths or periods by dividing vector length by the frequency
        # periods greater than half of the vector length are considered as 1 (no periodicity)
        periods = [0] * n
        ii = 0
        for i in ind:
            period = float(len(count_vector) + 1) / (i + 1)
            if period > len(count_vector) // 2:
                period = 1
            periods[ii] = period
            ii += 1

        # also return the whole thing as Series
        d = [float(len(count_vector) + 1) / (i + 1) for i in range(len(fft_abs))]
        fft_signal_series = pd.Series(data=fft_abs, index=d)

        return fft_signal_series, periods, signals, scales
