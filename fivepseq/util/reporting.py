"""
Classes and functions in this module will take care of fivepseq output: directory handling, writing of reports, etc.
"""
import fivepseq.config

class FivepseqOut:
    def __init__(self, output_dir):
        """

        :param output_dir:
        """

    def write_counts(self, count_array, file_name):
        """
        Writes the given count_array to the given file
        :param count_array:
        :param file_name:
        :return:
        """
        try:
            f = open(fivepseq.config.out_dir, "w")
        except Exception as e:
            error_message = "Problem opening the file %s for writing. Reason: %s" \
                            % (file_name, e.message)
            fivepseq.config.logger.error(error_message)
            raise Exception(error_message)
        for count in count_array:
            f.write(count)
        f.close()


