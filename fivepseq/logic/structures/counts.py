from math import floor

import plastid
from preconditions import preconditions

from fivepseq import config
from logic.structures.alignment import Alignment
from logic.structures.annotation import Annotation


class FivePSeqCounts:
    START = "start"
    TERM = "termination"
    FULL_LENGTH = "full_length"

    alignment = None
    annotation = None

    # @preconditions(lambda alignment: isinstance(alignment, Alignment),
    #                lambda annotation: isinstance(annotation, Annotation))
    def __init__(self, alignment, annotation):
        """
        Initializes a FivePSeqCounts object with Alignment and Annotation instances.
        :param alignment: Alignment type object
        :param annotation: Annotation type object
        """
        self.alignment = alignment
        self.annotation = annotation
        config.logger.debug("Initiated a FivePSeqCounts object with aligment from file %s and annotation from file %s."
                            % (alignment.alignment_file.filename, annotation.file_path))

    def get_counts(self, span_size, region=TERM):
        """
        Returns arrays of read counts within regions spanning the given position, with a span_size set in the annotation file.
        :param span_size: int: Specifies how many nucleotides to span around the specified region
        :param region: str: Specifies the region of the transcript to span for count vector generation
        :return: [[int]]: array of counts arrays of 5' mapping counts per position of the specified region of each transcript
        """

        config.logger.info("Retrieving counts from transcripts (span size :%d, position: %s)...."
                           % (span_size, region))

        count = self.annotation.transcript_count
        count_vectors = [None] * count
        progress_bar = ""
        i = 1
        transcript_generator = self.annotation.yield_transcripts(span_size)

        for transcript in transcript_generator:
            if i % 100 == 0:
                progress_bar += "#"
                print "\r>>Transcript count: %d\t%s\t%d%s" % (i, progress_bar, floor(100*(i-1)/count), '%'),
            try:
                count_vector = self.get_transcript_counts(transcript, span_size)
            except Exception as e:
                raise e
            count_vectors[i-1] = count_vector
            i += 1
        config.logger.info("Finished retrieving counts")
        return count_vectors

    @preconditions (lambda transcript: isinstance(transcript, plastid.genomics.roitools.Transcript),
                    lambda span_size: isinstance(span_size, int),
                    lambda region: isinstance(region, str))
    def get_transcript_counts(self, transcript, span_size, region = FULL_LENGTH):
        """
        Returns the array of counts per each transcript position within the given spanning region.
        :param region: str: Specifies the region of the transcript to span for count vector generation
        :param transcript: plastid.Transcript: The transcript to return the counts for: is is already spanned with the specified span_size
        :param span_size: int: Specifies how many nucleotides to span around the specified region
        :return: [int]: array of 5' mapping counts per position of the specified transcript region
        """
        try:
            count_vector = transcript.get_counts(self.alignment.bam_array)
            if region == self.FULL_LENGTH:
                return count_vector
            elif region == self.START:
                count_vector = count_vector[:2 * span_size]
            elif region == self.TERM:
                count_vector = count_vector[-(2 * span_size):]
            else:
                error_message = "Problem retrieving the count vector for the transcript %s. " \
                                "Invalid region \"%s\" specified: should be one of (%s, %s, %s)." \
                                % (transcript.get_name(), region, self.FULL_LENGTH, self.START, self.TERM)
                config.logger.error(error_message)
                raise ValueError(error_message)
        except Exception as e:
            error_message = "Problem retrieving the count vector for the transcript %s. Reason:%s" % (
                transcript.get_name(), e.message)
            config.logger.error(error_message)
            raise Exception(error_message)

        return count_vector
