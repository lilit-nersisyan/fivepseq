"""
Classes and functions in this module will take care of fivepseq output: directory handling, writing of reports, etc.
"""
import os
import numpy as np
from fivepseq import config
from preconditions import preconditions

import fivepseq.config


class FivepseqOut:
    output_dir = None

    def __init__(self, output_dir):
        """

        :param output_dir:
        """
        self.output_dir = output_dir



    @preconditions(lambda count_array: isinstance(count_array, list),
               lambda file_name: isinstance(file_name, str))
    def write_counts(self, count_array, file_name):
        """
        Writes the given count_array to the given file
        :param count_array:
        :param file_name:
        :return:
        """
        try:
            f = open(os.path.join(self.output_dir, file_name), "w")
        except Exception as e:
            error_message = "Problem opening the file %s for writing. Reason: %s" \
                            % (file_name, e.message)
            config.logger.error(error_message)
            raise Exception(error_message)
        config.logger.info("Storing count array in file %s..." % f)
        for counts in count_array:
            f.write('\t'.join(map(str, counts)) + '\n')
        f.close()
        config.logger.info("File %s saved" % f)

