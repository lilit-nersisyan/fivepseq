"""
Classes and functions in this module will take care of fivepseq output: directory handling, writing of reports, etc.
"""
import os
from pandas import DataFrame

from fivepseq import config
from preconditions import preconditions


class FivePSeqOut:
    output_dir = None

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

    @preconditions(lambda df: isinstance(df, DataFrame),
                   lambda file_name: isinstance(file_name, str))
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
                file.write('\t'.join(["ID", "gene", "cds_start", "cds_end", "type"]) + "\n")
                for transcript in transcript_assembly:
                    file.write('\t'.join([transcript.attr["ID"],
                                          transcript.attr["gene_id"],
                                          str(transcript.attr["cds_genome_start"]),
                                          str(transcript.attr["cds_genome_end"]),
                                          transcript.attr["type"]]) + '\n')
        except Exception as e:
            error_message = "Problem opening the file %s for writing. Reason: %s" \
                            % (file_name, e.message)
            config.logger.error(error_message)
            raise Exception(error_message)

