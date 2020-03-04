
from scipy import stats, special
import tasklogger
import pandas as pd
import numpy as np
import os
import sys

from . import utils


class SCYN:

    def __init__(self, seq='single-end', bin_len=500, ref='hg19', reg='*.bam', mapq=40, verbose=1):
        self.seq = seq
        self.bin_len = bin_len
        self.ref = ref
        self.reg = reg
        self.mapq = mapq
        self.verbose = verbose
        self._check_params()
        self.cnv = None
        self.meta_info = None
        self.segments = None
        self.bin_info = None
        tasklogger.set_level(verbose)


    def _check_params(self):
        """Check SCYN parameters

        This allows us to fail early - otherwise certain unacceptable
        parameter choices, such as n='10.5', would only fail after
        minutes of runtime.

        Raises
        ------
        ValueError : unacceptable choice of parameters
        """
        utils.check_in(['single-end', 'paired-end'], seq=self.seq)
        utils.check_int(bin_len=self.bin_len)
        utils.check_positive(bin_len=self.bin_len)
        utils.check_in(['hg19', 'hg38'], ref=self.ref)
        utils.check_int(mapq=self.mapq)
        utils.check_positive(mapq=self.mapq)


    def set_params(self, **params):
        """Set the parameters of SCYN.

        Any parameters not given as named arguments will be left at their
        current value.

        Parameters
        ----------

        seq : string, optional, default: single-end
            The reads type: single-end or paired-end

        bin_len : int, optional, default: 500
            The bin length, default is 500K
        
        ref : string, optional, default: hg19
            The reference genome version: hg19 or hg38

        reg : string, optional, default: *.bam
            The regular expression to match all BAM files in your input directory.
            For example, "*.bam" will match all BAM files ended with '.bam'
        
        mapq : int, optional, default: 40
            The mapping quality cutoff when calculating the reads coverage
        
        K : int, optional, default: 10
            The predifined changepoints number for all chromosomes

        verbose : `int` or `boolean`, optional, default: 1
            If `True` or `> 0`, print log messages

        Returns
        -------
        self
        """

        # parameters
        if 'seq' in params and params['seq'] != self.seq:
            self.seq = params['seq']
            del params['seq']
        if 'bin_len' in params and params['bin_len'] != self.bin_len:
            self.bin_len = params['bin_len']
            del params['bin_len']
        if 'ref' in params and params['ref'] != self.ref:
            self.ref = params['ref']
            del params['ref']
        if 'reg' in params and params['reg'] != self.reg:
            self.reg = params['reg']
            del params['reg']
        if 'mapq' in params and params['mapq'] != self.mapq:
            self.mapq = params['mapq']
            del params['mapq']
        if 'verbose' in params:
            self.verbose = params['verbose']
            tasklogger.set_level(self.verbose)
            del params['verbose']

        self._check_params()
        return self

    def _cal_cnv(self, ref, Y, Ynor):
        all_break_points = dict.fromkeys(np.unique(ref['seqnames']))
        all_cnv = pd.DataFrame(columns=Y.columns)
        keys = list(all_break_points.keys())
        if 'chr' in keys[0]:
            keys = sorted(list(keys), key=lambda x: int(x[3:]))
        for chrom in keys:
            indexes = np.where(ref['seqnames'] == chrom)
            # tasklogger.log_info('Calculate CNV for chromosome ' + chrom+'...')
            breaks, cnv = self._cal_cnv_for_each_chrom(Y.iloc[indexes], Ynor.iloc[indexes])
            all_break_points[chrom] = breaks
            all_cnv = all_cnv.append(cnv)
        self.cnv = all_cnv
        self.segments = pd.DataFrame(columns=['Chromosome', 'Start_bin_no', 'End_bin_no'])
        for chrom in keys:
            break_points = all_break_points[chrom]
            for i in range(len(break_points) - 1):
                start = break_points[i]
                end = break_points[i+1] -1
                self.segments = self.segments.append([{
                    'Chromosome': chrom,
                    'Start_bin_no': start,
                    'End_bin_no': end
                }], ignore_index=True)
        return self

    def _get_segments(self, break_points):
        segments = np.zeros((len(break_points) -1, 2))
        for i in range(len(break_points) - 1):
            start = break_points[i]
            end = break_points[i+1] - 1
            segments[i] = [start, end]
        return segments

    def _cal_mbic(self, Y, nor_Y, segs):
        nor_Y.columns = Y.columns
        nor_Y.index = Y.index
        Y0 = Y.copy()
        nor_Y0 = nor_Y.copy()
        index_Y = np.column_stack(np.where(Y0 <= 20))
        for index in index_Y:
            Y0.iat[index[0], index[1]] = 20
        index_nor_Y = np.column_stack(np.where(nor_Y0 <= 20))
        for index in index_nor_Y:
            nor_Y0.iat[index[0], index[1]] = 20
        bin_num = Y.shape[0]
        sample_num = Y.shape[1]

        avg_mat = np.zeros((segs.shape[0], sample_num))
        rate_mat = np.zeros((segs.shape[0], sample_num))
        sum_Y_mat = np.zeros((segs.shape[0], sample_num))
        sum_nor_Y_mat = np.zeros((segs.shape[0], sample_num))
        sum_Y_mat0 = np.zeros((segs.shape[0], sample_num))
        sum_nor_Y_mat0 = np.zeros((segs.shape[0], sample_num))
        cumsum_Y = np.cumsum(Y).values
        cumsum_nor_Y = np.cumsum(nor_Y).values
        cumsum_Y0 = np.cumsum(Y0).values
        cumsum_nor_Y0 = np.cumsum(nor_Y0).values
        for index in range(segs.shape[0]):
            row = segs[index]
            start = int(row[0])
            end = int(row[1])
            if start == 0:
                sum_Y_mat[index] = cumsum_Y[end]
                sum_nor_Y_mat[index] = cumsum_nor_Y[end]
                sum_Y_mat0[index] = cumsum_Y0[end]
                sum_nor_Y_mat0[index] = cumsum_nor_Y0[end]
            else:
                sum_Y_mat[index] = cumsum_Y[end] - cumsum_Y[start - 1]
                sum_nor_Y_mat[index] = cumsum_nor_Y[end] - cumsum_nor_Y[start - 1]
                sum_Y_mat0[index] = cumsum_Y0[end] - cumsum_Y0[start - 1]
                sum_nor_Y_mat0[index] = cumsum_nor_Y0[end] - cumsum_nor_Y0[start - 1]
            avg_mat[index] = np.round(sum_Y_mat0[index]/(end-start+1))
            # print(end-start+1)
            rate_mat[index] = np.round(sum_Y_mat0[index]/sum_nor_Y_mat0[index] * 2)
        k1 = 3/2
        k2 = 2.27
        i = segs.shape[0] - 1
        # print(avg_mat.shape, rate_mat.shape, sum_Y_mat.shape, sum_nor_Y_mat.shape)
        rate_carriers = rate_mat[1:] - rate_mat[:i]
        rate_carriers[rate_carriers != 0] = 1
        avg_carriers = avg_mat[1:] - avg_mat[:i]
        total_carriers = avg_carriers * rate_carriers
        max_sum = np.max(
            np.sum(np.power(total_carriers, 2), axis=1))
        M = np.sum(rate_carriers)
        pi = M / (sample_num * (i))
        loglikeij = 0.0

        for row_index in np.arange(sum_Y_mat.shape[0]):
            row_sum_Y = sum_Y_mat[row_index]
            row_sum_nor_Y = sum_nor_Y_mat[row_index]
            row_sum_Y[row_sum_nor_Y < 20] = 20
            row_sum_nor_Y[row_sum_nor_Y < 20] = 20
            loglikeij += np.sum((1 - np.round(row_sum_Y / row_sum_nor_Y * 2) / 2) * row_sum_nor_Y
                                + np.log((np.round(row_sum_Y / row_sum_nor_Y * 2) + 1e-04) / 2.0001) * row_sum_Y)
        term1 = loglikeij
        # print(sys.argv[3], term1)
        if M == 0:
            term2 = 0
        else:
            term2 = -M / 2 * np.log(2 * loglikeij/M)
        term3 = -np.log(special.binom(bin_num, i))
        term4 = -M / 2
        if M == 0 or max_sum == 0:
            term5 = 0
        else:
            term5 = -np.sum(np.log(max_sum))
        term6 = -(i) * (k1 - k2)
        # print(pi)
        if pi == 0 or pi == 1:
            term7 = 0
        else:
            term7 = (M * np.log(pi) + (sample_num * (i) - M)
                    * np.log(1 - pi))
        mbic = term1 + term2 + term3 + term4 + term5 + term6 + term7
        return mbic
    
    def _get_break_points_for_each_k(self, paths, i):
        break_points = []
        break_points.append(paths.shape[1])

        j = paths.shape[1] - 1
        while i >= 0:
            break_points.append(paths[i][j])
            j = paths[i][j] - 1
            i = i - 1
        if break_points[-1] != 0:
            break_points.append(0)
        break_points.reverse()
        return break_points

    def _cal_cnv_for_each_chrom(self, Y, nor_Y):
        nor_Y.columns = Y.columns
        nor_Y.index = Y.index
        bin_num = Y.shape[0]
        # mBIC = np.zeros((bin_num - 1, bin_num))
        # paths = np.zeros((bin_num - 1, bin_num), dtype=np.int64)
        mBIC = np.zeros((30, bin_num))
        paths = np.zeros((30, bin_num), dtype=np.int64)
        loglikeijs = np.zeros((bin_num, bin_num))
        cumsum_Y = np.cumsum(Y).values
        cumsum_nor_Y = np.cumsum(nor_Y).values
        
        
        for i in range(bin_num):
            for j in range(i, bin_num):
                if i == 0:
                    row_sum_Y = cumsum_Y[j]
                    row_sum_nor_Y = cumsum_nor_Y[j]
                else:
                    row_sum_Y = cumsum_Y[j] - cumsum_Y[i - 1]
                    row_sum_nor_Y = cumsum_nor_Y[j] - cumsum_nor_Y[i - 1]
                row_sum_Y[row_sum_nor_Y < 20] = 20
                row_sum_nor_Y[row_sum_nor_Y < 20] = 20
                loglikeij = np.sum((1 - np.round(row_sum_Y / row_sum_nor_Y * 2) / 2) * row_sum_nor_Y
                                    + np.log((np.round(row_sum_Y / row_sum_nor_Y * 2) + 1e-04) / 2.0001) * row_sum_Y)
                loglikeijs[i][j] = loglikeij
        k1 = 3/2
        k2 = 2.27
        i = 0
        count = 0
        while True:
            # if i >= bin_num - 1:
            if i >= 30 or i >= bin_num - 1:
                break
            for j in range(i + 1, bin_num):
                max_mbic = None

                # cal all possible mbics for current bin size
                for index in np.arange(i, j):
                    if i == 0:
                        term1 = loglikeijs[i][index] + loglikeijs[index + 1][j]
                    else:
                        term1 = mBIC[i-1][index] + loglikeijs[index + 1][j]

                    term2 = -np.log(special.binom(j+1, i + 1))
                    term3 = -(i + 1) * (k1 - k2)
                    mbic = term1 + term2 + term3
                    if index == i or mbic >= max_mbic:
                        mBIC[i][j] = mbic
                        paths[i][j] = index + 1
                        max_mbic = mbic
            last_col = mBIC[:, -1]
            if i != 0 and count == 0:
                if last_col[i] < last_col[i - 1]:
                    count = count + 1
            if count != 0:
                count = count + 1
            if count > 5:
                break
            i = i + 1
        
        max_k = i
        break_points_of_all_k = []
        mbic_of_all_k = []
        for i in range(max_k):
            break_points = self._get_break_points_for_each_k(paths, i)
            segs = self._get_segments(break_points)
            mbic = self._cal_mbic(Y, nor_Y, segs)
            break_points_of_all_k.append(break_points)
            mbic_of_all_k.append(mbic)
        # print(max(mbic_of_all_k))
        max_mbic_index = mbic_of_all_k.index(max(mbic_of_all_k))
        # print(max_mbic_index)
        break_points = break_points_of_all_k[max_mbic_index]
        last_col = mBIC[:, -1]
        print(last_col[max_mbic_index])
        # cal cnv for each segment
        index_Y = np.column_stack(np.where(Y <= 20))
        for index in index_Y:
            Y.iat[index[0], index[1]] = 20
        index_nor_Y = np.column_stack(np.where(nor_Y <= 20))
        for index in index_nor_Y:
            nor_Y.iat[index[0], index[1]] = 20
        cumsum_Y = np.cumsum(Y).values
        cumsum_nor_Y = np.cumsum(nor_Y).values
        cnv = np.zeros(Y.shape, dtype=np.int64)
        for i in range(len(break_points) - 1):
            start = break_points[i]
            end = break_points[i+1] - 1
            if start == 0:
                sum_Y = cumsum_Y[end]
                sum_nor_Y = cumsum_nor_Y[end]
            else:
                sum_Y = cumsum_Y[end] - cumsum_Y[start - 1]
                sum_nor_Y = cumsum_nor_Y[end] - cumsum_nor_Y[start - 1]
            rate = np.round((sum_Y / sum_nor_Y) * 2)
            rate[rate > 14] = 14
            
            for j in np.arange(start, end+1):
                cnv[j] = rate
        cnv = pd.DataFrame(cnv, columns=Y.columns, index=Y.index)
        return break_points, cnv


    def call(self, bam_dir, out_dir):
        """call CNV for each chromosome

        Parameters
        ----------
        bam_dir : directory path which contains all BAM files

        out_dir : the output directory
        Returns
        -------
        self 
        """
        Y_path = os.path.join(out_dir, 'temp.Y.csv')
        nor_Y_path = os.path.join(out_dir, 'temp.norY.csv')
        ref_path = os.path.join(out_dir, 'temp.ref.csv')
        gini_path = os.path.join(out_dir, 'temp.gini.csv')
        ploidy_path = os.path.join(out_dir, 'temp.ploidy.csv')
        scope_path = os.path.join(out_dir, 'run-scope.R')
        utils.write_scope(scope_path)
        command = 'Rscript {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}'.format(scope_path, bam_dir, Y_path, ref_path, gini_path, ploidy_path, nor_Y_path,self.seq, self.reg, self.ref, self.mapq, self.bin_len)
        code = os.system(command)
        if code != 0:
            sys.exit(1)
        tasklogger.log_start('SCYN')
        Y = pd.read_csv(Y_path, index_col=0)
        nor_Y = pd.read_csv(nor_Y_path, index_col=0)
        ref = pd.read_csv(ref_path, index_col=0)
        gini = pd.read_csv(gini_path, index_col=0)
        ploidy = pd.read_csv(ploidy_path, index_col=0)
        self.meta_info = pd.DataFrame(index=['c_gini', 'c_ploidy'], columns=Y.columns)
        self.meta_info.loc['c_gini'] = gini.T.iloc[0].values
        self.meta_info.loc['c_ploidy'] = ploidy.T.iloc[0].values
        self.meta_info = self.meta_info.T
        self._cal_cnv(ref, Y, nor_Y)
        self.bin_info = ref

        # clean up temp files
        utils.clean_up([Y_path, nor_Y_path, ref_path,
                        gini_path, ploidy_path, scope_path])

        tasklogger.log_complete('SCYN')

        return self

        

    
