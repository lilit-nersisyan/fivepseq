"""
Classes and functions in this module will take care of fivepseq output: directory handling, writing of reports, etc.
"""
import os
import pandas as pd
from pandas import DataFrame

from fivepseq import config
from preconditions import preconditions


class FivePSeqOut:
    output_dir = None
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

    def __init__(self, output_dir):
        """
        Initialized a reporter class object with the given directory.
        All further requests to store data will lead to files saved in the output directory defined in self.

        :param output_dir: str: the path to the fivepseq output directory
        """
        self.output_dir = output_dir

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
        try:
            f = open(os.path.join(self.output_dir, file_name), "w")
        except Exception as e:
            error_message = "Problem opening the file %s for writing. Reason: %s" \
                            % (file_name, e.message)
            config.logger.error(error_message)
            raise Exception(error_message)
        config.logger.info("Storing vector list to file %s..." % f)
        for counts in vector_list:
            f.write('\t'.join(map(str, counts)) + '\n')
        f.close()
        config.logger.info("File %s saved" % f)


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
        try:
            f = open(os.path.join(self.output_dir, file_name), "w")
        except Exception as e:
            error_message = "Problem opening the file %s for writing. Reason: %s" \
                            % (file_name, e.message)
            config.logger.error(error_message)
            raise Exception(error_message)
        config.logger.info("Storing vector list to file %s..." % f)
        for key  in vector_dict.keys():
            counts = vector_dict.get(key)
            f.write(key + "\t")
            f.write('\t'.join(map(str, counts)) + '\n')
        f.close()
        config.logger.info("File %s saved" % f)

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
        try:
            f = open(os.path.join(self.output_dir, file_name), "w")
        except Exception as e:
            error_message = "Problem opening the file %s for writing. Reason: %s" \
                            % (file_name, e.message)
            config.logger.error(error_message)
            raise Exception(error_message)
        config.logger.info("Storing vector list to file %s..." % f)
        for key in my_dict.keys():
            count = my_dict.get(key)
            f.write(key + "\t")
            f.write(str(count) + '\n')
        f.close()
        config.logger.info("File %s saved" % f)

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
        try:
            f = open(os.path.join(self.output_dir, file_name), "w")
        except Exception as e:
            error_message = "Problem opening the file %s for writing. Reason: %s" \
                            % (file_name, e.message)
            config.logger.error(error_message)
            raise Exception(error_message)
        config.logger.info("Storing vector list to file %s..." % f)
        for i in my_vector:
            f.write(str(i) + "\n")
        f.close()
        config.logger.info("File %s saved" % f)

    #@preconditions(lambda df: isinstance(df, pd.core.frame.DataFrame),
    #               lambda file_name: isinstance(file_name, str))
    def write_df_to_file(self, df, file_name):
        """
        Stores the given dataframe as tab-delimited text file
        The file is stored in the output directory defined in self.

        :param df: pandas.DataFrame
        :param file_name: str: the name of the file (without directory path)

        :return:
        """
        try:
            file = os.path.join(self.output_dir, file_name)
        except Exception as e:
            error_message = "Problem opening the file %s for writing. Reason: %s" \
                            % (file_name, e.message)
            config.logger.error(error_message)
            raise Exception(error_message)
        df.to_csv(file, sep="\t", header=True)

        config.logger.info("Dataframe saved to file %s " % file)

    def write_series_to_file(self, data_series, file_name):
        """
        Stores the indices and values of the given data series as two tab-delimited columns.
        The file is stored in the output directory defined in self.

        :param data_series: pandas.Series
        :param file_name: str: the name of the file without directory path

        :return: nothing to return
        """
        try:
            file = os.path.join(self.output_dir, file_name)
        except Exception as e:
            error_message = "Problem opening the file %s for writing. Reason: %s" \
                            % (file_name, e.message)
            config.logger.error(error_message)
            raise Exception(error_message)
        data_series.to_csv(file, sep="\t")

    def write_transcript_assembly_to_file(self, transcript_assembly, file_name):

        try:
            with open(os.path.join(self.output_dir, file_name), 'w') as file:
                # print the transcript attributes to the file
                file.write('\t'.join(["ID", "gene", "chr", "cds_start", "cds_end", "type"]) + "\n")
                for transcript in transcript_assembly:
                    file.write('\t'.join([transcript.attr["ID"],
                                          transcript.attr["gene_id"],
                                          str(transcript.chrom),
                                          str(transcript.attr["cds_genome_start"]),
                                          str(transcript.attr["cds_genome_end"]),
                                          transcript.attr["type"]]) + '\n')
        except Exception as e:
            error_message = "Problem opening the file %s for writing. Reason: %s" \
                            % (file_name, e.message)
            config.logger.error(error_message)
            raise Exception(error_message)

    def get_file_path(self, file_name):
        return os.path.join(self.output_dir, file_name)