from math import floor

from preconditions import preconditions

from fivepseq import config
from logic.structures.alignment import Alignment
from logic.structures.annotation import Annotation


class FivePSeqCounts:
    START = "start"
    TERM = "termination"

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

    def get_counts(self, span_size, position=TERM):
        """
        Returns arrays of read counts within regions spanning the given position, with a span_size set in the annotation file.
        :param span_size:
        :param position: either start or end of the transcript, or a specific position (later for amino acids)
        :return:
        """

        config.logger.info("Retrieving counts from transcripts (span size :%d, position: %s)...."
                           % (span_size, position))

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

    def get_transcript_counts(self, transcript, span_size):
        try:
            count_vector = transcript.get_counts(self.alignment.bam_array)
            if len(count_vector) >= 2 * span_size:
                count_vector = count_vector[:2 * span_size]
        except Exception as e:
            error_message = "Problem retrieving the count vector for the transcript %s. Reason:%s" % (
                transcript.get_name(), e.message)
            config.logger.error(error_message)
            raise Exception(error_message)

        return count_vector
