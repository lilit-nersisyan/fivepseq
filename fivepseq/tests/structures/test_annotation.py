import unittest

import plastid

from fivepseq.logic.structures.annotation import Annotation
from fivepseq.util.readers import AnnotationReader


class AnnotationTest(unittest.TestCase):
    ANNOT_FILE_PATH_VALID_F = "/proj/sllstore2017018/lilit/5pseq_microbiome/gff/Bacillus_subtilis.ASM69118v1.37.gff3"
    annotation = None
    def setUp(self):
        annotation_reader = AnnotationReader(self.ANNOT_FILE_PATH_VALID_F)
        self.annotation = annotation_reader.annotation

    def test_annotation_init(self):
        with self.assertRaises(ValueError):
            Annotation(None)

        self.assertIsNotNone(self.annotation.transcript_assembly)
        self.assertIsInstance(self.annotation.transcript_assembly,
                              plastid.readers.gff.GFF3_TranscriptAssembler)

