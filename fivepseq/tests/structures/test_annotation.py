import unittest

import numpy

from fivepseq.util.readers import AnnotationReader


class AnnotationTest(unittest.TestCase):
    ANNOT_FILE_PATH_VALID_F = "/proj/sllstore2017018/lilit/5pseq_microbiome/gff/Bacillus_subtilis.ASM69118v1.37.gff3"
    annotation = None

    def setUp(self):
        annotation_reader = AnnotationReader(self.ANNOT_FILE_PATH_VALID_F, break_after=10)
        self.annotation = annotation_reader.annotation

    def test_annotation_init(self):
        self.assertIsNotNone(self.annotation.transcript_count)

    def test_yield_transcripts(self):

        t0_array = []
        t1_array = []
        t2_array = []
        t2_array2 = []
        for t in self.annotation.yield_transcripts(1):
            t0_array.append(t)
        for t in self.annotation.yield_transcripts(10):
            t1_array.append(t)
        for t in self.annotation.yield_transcripts(20):
            t2_array.append(t)
        for t in self.annotation.yield_transcripts(20):
            t2_array2.append(t)
        self.assertEqual(t0_array[0].get_name(), t1_array[0].get_name())
        self.assertTrue(t0_array[0].spanning_segment.start == 0 or t0_array[0].spanning_segment.start > t1_array[
            0].spanning_segment.start)
        self.assertTrue(t1_array[0].spanning_segment.end + 10 == t2_array[0].spanning_segment.end)
        numpy.testing.assert_equal(t1_array, t1_array)
        numpy.testing.assert_equal(t2_array, t2_array2)
