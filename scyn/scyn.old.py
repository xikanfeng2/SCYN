
from scipy import stats, special
import tasklogger
import pandas as pd
import numpy as np
import os
import sys

from . import utils


class SCYN:

    def __init__(self, seq='single-end', bin_len=500, ref='hg19', reg='*.bam', mapq=40, K=10, verbose=1):
        self.seq = seq
        self.bin_len = bin_len
        self.ref = ref
        self.reg = reg
        self.mapq = mapq
        self.K = K
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
        utils.check_int(K=self.K)
        utils.check_positive(K=self.K)


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
        if 'K' in params and params['K'] != self.K:
            self.K = params['K']
            del params['K']
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
            tasklogger.log_info('Calculate CNV for chromosome ' + chrom+'...')
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
    
    def _cal_cnv_for_each_chrom(self, Y, nor_Y):
        nor_Y.columns = Y.columns
        nor_Y.index = Y.index
        index_Y = np.column_stack(np.where(Y <= 20))
        for index in index_Y:
            Y.iat[index[0], index[1]] = 20
        index_nor_Y = np.column_stack(np.where(nor_Y <= 20))
        for index in index_nor_Y:
            nor_Y.iat[index[0], index[1]] = 20
        bin_num = Y.shape[0]
        sample_num = Y.shape[1]
        mBIC = pd.DataFrame(columns=np.arange(bin_num), index=np.arange(self.K), dtype=float)
        paths = pd.DataFrame(columns=np.arange(bin_num),
                            index=np.arange(self.K))
        avg_mats = np.empty((self.K + 1, bin_num), dtype=object)
        rate_mats = np.empty((self.K + 1, bin_num), dtype=object)
        sum_Y_mats = np.empty((self.K+1, bin_num), dtype=object)
        sum_nor_Y_mats = np.empty((self.K+1, bin_num), dtype=object)
        cumsum_Y = np.cumsum(Y).values
        cumsum_nor_Y = np.cumsum(nor_Y).values
        init_avg_result = np.round(
            cumsum_nor_Y/np.arange(1, bin_num+1)[:, None])
        init_rate_result = np.round(cumsum_Y/cumsum_nor_Y * 2)
        
        for i in range(bin_num):
            avg_mats[0][i] = [init_avg_result[i]]
            rate_mats[0][i] = [init_rate_result[i]]
            sum_Y_mats[0][i] = [cumsum_Y[i]]
            sum_nor_Y_mats[0][i] = [cumsum_nor_Y[i]]

        k1 = 3/2
        k2 = 2.27
        
        for i in np.arange(self.K):
            print(i)
            for j in np.arange(bin_num):
                if j <= i:
                    continue
                max_mbic = 0.0

                # cal all possible mbics for current bin size
                for index in np.arange(i, j):
                    last_sum_row_Y = cumsum_Y[j] - cumsum_Y[index]
                    last_sum_row_nor_Y = cumsum_nor_Y[j] - cumsum_nor_Y[index]
                    last_avg_row = np.round(last_sum_row_nor_Y/(j-i))
                    last_rate_row = np.round(last_sum_row_Y/last_sum_row_nor_Y * 2)
                    previous_avg_mat = avg_mats[i, index]
                    avg_mat = np.append(previous_avg_mat, [last_avg_row], axis=0)
                    previous_rate_mat = rate_mats[i, index]
                    rate_mat = np.append(previous_rate_mat, [last_rate_row], axis=0)
                    previous_sum_Y_mat = sum_Y_mats[i, index]
                    sum_Y_mat = np.append(previous_sum_Y_mat, [last_sum_row_Y], axis=0)
                    previous_sum_nor_Y_mat = sum_nor_Y_mats[i, index]
                    sum_nor_Y_mat = np.append(previous_sum_nor_Y_mat, [last_sum_row_nor_Y], axis=0)

                    # cal mbic
                    rate_carriers = rate_mat[1:] - rate_mat[:i+1]
                    rate_carriers[rate_carriers != 0] = 1
                    avg_carriers = avg_mat[1:] - avg_mat[:i+1]
                    total_carriers = avg_carriers * rate_carriers
                    max_sum = np.max(
                        np.sum(np.power(total_carriers, 2), axis=1))
                    M = np.sum(rate_carriers)
                    pi = M / (sample_num * (i + 1))
                    loglikeij = 0.0

                    for row_index in np.arange(sum_Y_mat.shape[0]):
                        row_sum_Y = sum_Y_mat[row_index]
                        row_sum_nor_Y = sum_nor_Y_mat[row_index]
                        loglikeij += np.sum((1 - np.round(row_sum_Y / row_sum_nor_Y * 2) / 2) * row_sum_nor_Y
                                            + np.log((np.round(row_sum_Y / row_sum_nor_Y * 2) + 1e-04) / 2.0001) * row_sum_Y)
                    term1 = loglikeij
                    if M == 0 or term1 <= 0:
                        term2 = 0
                    else:
                        term2 = -M / 2 * np.log(2 * loglikeij/M)
                    term3 = -np.log(special.binom(j+1, i + 1))
                    term4 = -M / 2
                    if M == 0 or max_sum == 0:
                        term5 = 0
                    else:
                        term5 = -np.sum(np.log(max_sum))
                    term6 = -(i + 1) * (k1 - k2)
                    if pi == 0 or pi == 1:
                        term7 = 0
                    else:
                        term7 = (M * np.log(pi) + (sample_num * (i + 1) - M)
                                 * np.log(1 - pi))
                    mbic = term1 + term2 + term3 + term4 + term5 + term6 + term7
                    if np.isnan(mBIC.iloc[i, j]) or mbic >= max_mbic:
                        mBIC.iloc[i, j] = mbic
                        paths.iloc[i, j] = index + 1
                        sum_Y_mats[i+1, j] = sum_Y_mat
                        sum_nor_Y_mats[i+1, j] = sum_nor_Y_mat
                        avg_mats[i+1, j] = avg_mat
                        rate_mats[i+1, j] = rate_mat
                        max_mbic = mbic         
        
        # backtrack
        break_points = []
        break_points.append(bin_num)
        max_k = mBIC.iloc[:, -1].idxmax()
        i = max_k
        j = paths.shape[1] - 1
        while i >= 0:
            break_points.append(paths.values[i, j])
            j = paths.values[i, j] - 1
            i = i - 1
        if break_points[-1] != 0:
            break_points.append(0)
        break_points.reverse()
        print(break_points)
        # cal cnv for each segment
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
        tasklogger.log_complete('SCYN')

        return self

        

    
