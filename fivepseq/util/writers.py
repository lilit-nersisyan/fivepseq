"""
Classes and functions in this module will take care of fivepseq output: directory handling, writing of reports, etc.
"""
import logging
import os

from preconditions import preconditions

from fivepseq import config


class FivePSeqOut:
    output_dir = None
    conflict_mode = None
    # file names
    COUNT_TERM_FILE = "counts_TERM.txt"
    COUNT_START_FILE = "counts_START.txt"
    COUNT_FULL_FILE = "counts_FULL_LENGTH.txt"
    META_COUNT_TERM_FILE = "meta_counts_TERM.txt"
    META_COUNT_START_FILE = "meta_counts_START.txt"
    AMINO_ACID_PAUSES_FILE = "amino_acid_pauses.txt"
    LOCI_PAUSES_FILE = "loci_pauses.txt"
    TRANSCRIPT_ASSEMBLY_FILE = "transcript_assembly.txt"
    FRAME_COUNTS_TERM_FILE = "frame_counts_TERM.txt"
    FRAME_COUNTS_START_FILE = "frame_counts_START.txt"
    START_CODON_DICT_FILE = "start_codon_dict.txt"
    TERM_CODON_DICT_FILE = "term_codon_dict.txt"
    CANONICAL_TRANSCRIPT_INDEX_FILE = "canonical_transcript_index.txt"
    TRANSCRIPT_DESCRIPTORS_FILE = "transcript_descriptors.txt"
    DATA_SUMMARY_FILE = "data_summary.txt"
    FRAME_STATS_DF_FILE = "frame_stats.txt"
    FFT_STATS_DF_FILE = "fft_stats.txt"
    FFT_SIGNALS_START = "fft_signals_start.txt"
    FFT_SIGNALS_TERM = "fft_signals_term.txt"

    FAILED_COUNT_FILES_LIST = "failed_count_files_list.txt"
    BAM_SUCCESS_SUMMARY = "count_run_summary.txt"

    logger = logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER)
    count_success = None

    def __init__(self, output_dir, conflict_mode=config.ADD_FILES):
        """
        Initialized a reporter class object with the given directory.
        All further requests to store data will lead to files saved in the output directory defined in self.

        :param output_dir: str: the path to the fivepseq output directory
        """
        self.output_dir = output_dir
        self.conflict_mode = conflict_mode

    def open_file_for_writing(self, file_name):
        """
        Prepares the given file for writing in self output directory.
        If the specified file already exists, the conflict is handled according to self conflict handling mode.
        If the conflict handling mode equals:
            -config.ADD_FILES: the file is not written as it already exists.
            -config.OVERWRITE: The existing file is removed. If not possible an Error is raised.
            -config.ALT_DIR: Such cases should not occur, as this a brand new directory. If they occur, an Error should be raised.

        :param file_name: str: the name of the file to be stored in the output directory.
        :return: the file object opened for writing
        """

        path = os.path.join(self.output_dir, file_name)

        if os.path.exists(path):

            if self.conflict_mode == config.ADD_FILES:
                self.logger.info("Skipping existing file %s." % path)
                return None

            elif self.conflict_mode == config.ALT_DIR:
                raise Exception("The file %s already exists in the output directory %s.\n "
                                "Under the %s conflict handling mode this should be a brand new empty directory.\n"
                                "If you are here then you're writing to the same file twice"
                                % (file_name, self.output_dir, self.conflict_mode))
            else:
                try:
                    os.remove(path)
                    self.logger.debug("Removed existing file %s for overwriting" % path)
                except Exception as e:
                    self.logger.error("Could not remove existing file %s for overwriting. "
                                      "Reason: %s" % (path, str(e)))

        try:
            f = open(path, "w")
        except Exception as e:
            error_message = "Problem opening the file %s for writing. Reason: %s" \
                            % (file_name, e.message)
            self.logger.error(error_message)
            raise Exception(error_message)

        return f

    @preconditions(lambda vector_list: isinstance(vector_list, list),
                   lambda file_name: isinstance(file_name, str))
    def write_vector_list(self, vector_list, file_name):
        """
        Writes the given vector list to the given file.
        The file is stored in the output directory defined in self.

        :param vector_list: [[int]]: a list of int vectors
        :param file_name: str: the name of the file (without directory path)

        :return: nothing to return
        """
        f = self.open_file_for_writing(file_name)
        if f is not None:
            self.logger.info("Storing vector list to file %s..." % f)
            for counts in vector_list:
                f.write('\t'.join(map(str, counts)) + '\n')
            f.close()
            self.logger.debug("File %s saved" % f)

    @preconditions(lambda vector_dict: isinstance(vector_dict, dict),
                   lambda file_name: isinstance(file_name, str))
    def write_vector_dict(self, vector_dict, file_name):
        """
        Writes the given vector dictionary to the given file.
        The file is stored in the output directory defined in self.

        :param vector_dict: {[int]}: a dict of int vectors
        :param file_name: str: the name of the file (without directory path)

        :return: nothing to return
        """
        f = self.open_file_for_writing(file_name)
        if f is not None:
            self.logger.info("Storing vector list to file %s..." % f)
            for key in vector_dict.keys():
                counts = vector_dict.get(key)
                f.write(key + "\t")
                f.write('\t'.join(map(str, counts)) + '\n')
            f.close()
            self.logger.debug("File %s saved" % f)

    @preconditions(lambda my_dict: isinstance(my_dict, dict),
                   lambda file_name: isinstance(file_name, str))
    def write_dict(self, my_dict, file_name):
        """
        Writes the given dictionary to the given file.
        The file is stored in the output directory defined in self.

        :param my_dict: {int}: a dict of codon : count pairs
        :param file_name: str: the name of the file (without directory path)

        :return: nothing to return
        """
        f = self.open_file_for_writing(file_name)
        if f is not None:
            self.logger.info("Storing vector list to file %s..." % f)
            for key in my_dict.keys():
                count = my_dict.get(key)
                f.write(key + "\t")
                f.write(str(count) + '\n')
            f.close()
            self.logger.debug("File %s saved" % f)

    @preconditions(lambda my_vector: isinstance(my_vector, list),
                   lambda file_name: isinstance(file_name, str))
    def write_vector(self, my_vector, file_name):
        """
        Writes the given vector to the given file.
        The file is stored in the output directory defined in self.

        :param my_vector: [int]: a vector of counts/indices
        :param file_name: str: the name of the file (without directory path)

        :return: nothing to return
        """
        f = self.open_file_for_writing(file_name)
        if f is not None:
            self.logger.info("Storing vector list to file %s..." % f)
            for i in my_vector:
                f.write(str(i) + "\n")
            f.close()
            self.logger.debug("File %s saved" % f)

    # @preconditions(lambda df: isinstance(df, pd.core.frame.DataFrame),
    #               lambda file_name: isinstance(file_name, str))
    def write_df_to_file(self, df, file_name):
        """
        Stores the given dataframe as tab-delimited text file
        The file is stored in the output directory defined in self.

        :param df: pandas.DataFrame
        :param file_name: str: the name of the file (without directory path)

        :return:
        """
        f = self.open_file_for_writing(file_name)
        if f is not None:
            df.to_csv(f, sep="\t", header=True)

            self.logger.debug("Dataframe saved to file %s " % f)

    def write_series_to_file(self, data_series, file_name):
        """
        Stores the indices and values of the given data series as two tab-delimited columns.
        The file is stored in the output directory defined in self.

        :param data_series: pandas.Series
        :param file_name: str: the name of the file without directory path

        :return: nothing to return
        """
        f = self.open_file_for_writing(file_name)
        if f is not None:
            data_series.to_csv(f, sep="\t")

    def write_transcript_assembly_to_file(self, transcript_assembly, file_name):
        f = self.open_file_for_writing(file_name)
        if f is not None:
            f.write('\t'.join(["ID", "gene", "chr", "cds_start", "cds_end", "type"]) + "\n")
            for transcript in transcript_assembly:
                f.write('\t'.join([self.get_transcript_attr(transcript, "ID"),
                                   self.get_transcript_attr(transcript, "Name"),
                                   str(transcript.chrom),
                                   str(self.get_transcript_attr(transcript, "cds_genome_start")),
                                   str(self.get_transcript_attr(transcript, "cds_genome_end")),
                                   self.get_transcript_attr(transcript, "type")]) + '\n')

            f.close()

    @staticmethod
    def get_transcript_attr(transcript, key):
        if key in transcript.attr.keys():
            return transcript.attr[key]
        return ""

    def get_file_path(self, file_name):
        return os.path.join(self.output_dir, file_name)

    def sanity_check_for_counts(self):
        """"
        Checks if all the count files are present in the output directory.
        Stores absent files in a text file.

        :return: returns True if all the files are in place and False otherwise
        """
        if self.count_success is not None:
            return self.count_success

        if os.path.exists(self.get_file_path(self.FAILED_COUNT_FILES_LIST)):
            try:
                os.remove(self.get_file_path(self.FAILED_COUNT_FILES_LIST))
            except Exception as e:
                logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).error\
                    ("Problem removing existing failed files report: %s" % str(e))

        failed_files = []

        necessary_files = [self.COUNT_TERM_FILE,
                           self.COUNT_START_FILE,
                           self.COUNT_FULL_FILE,
                           self.META_COUNT_TERM_FILE,
                           self.META_COUNT_START_FILE,
                           self.AMINO_ACID_PAUSES_FILE,
                           self.TRANSCRIPT_ASSEMBLY_FILE,
                           self.FRAME_COUNTS_TERM_FILE,
                           self.FRAME_COUNTS_START_FILE,
                           self.CANONICAL_TRANSCRIPT_INDEX_FILE,
                           self.TRANSCRIPT_DESCRIPTORS_FILE,
                           self.DATA_SUMMARY_FILE,
                           self.FRAME_STATS_DF_FILE,
                           self.FFT_STATS_DF_FILE,
                           self.FFT_SIGNALS_START,
                           self.FFT_SIGNALS_TERM]

        for f in necessary_files:
            if not os.path.exists(self.get_file_path(f)):
                failed_files.append(f)

        if len(failed_files) > 0:
            self.write_vector(failed_files, self.FAILED_COUNT_FILES_LIST)
            return False

        return True
