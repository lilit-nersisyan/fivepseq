import sys
import unittest
from StringIO import StringIO
from contextlib import contextmanager
from fivepseq.util.readers import *
from fivepseq.util import readers

from fivepseq.tests import test_files


class FileCheckTest(unittest.TestCase):

    def test_check_file_validity_valid_f(self):
        check_file_path(test_files.ALIGNMENT_FILE_PATH_VALID_F)

    def test_check_file_validity_d(self):
        with self.assertRaises(IOError):
            check_file_path(test_files.ALIGNMENT_FILE_PATH_D)

    def test_check_file_validity_invalid_f(self):
        with self.assertRaises(IOError):
            check_file_path(test_files.ALIGNMENT_FILE_PATH_INVALID_F)

    def test_get_file_compression(self):
        self.assertEqual(readers.get_file_compression("test.gz"), COMPRESSION_GZ)
        self.assertEqual(readers.get_file_compression("test.bz"), COMPRESSION_BZ)
        self.assertEqual(readers.get_file_compression("test.bz.aaa"), COMPRESSION_None)
        self.assertEqual(readers.get_file_compression("test.aaa"), COMPRESSION_None)
        self.assertEqual(readers.get_file_compression("test"), COMPRESSION_None)

    def test_get_file_extension(self):
        self.assertEqual(readers.get_file_extension("test.ext1"), "ext1")
        self.assertEqual(readers.get_file_extension("test.ext1.gz"), "ext1")
        self.assertEqual(readers.get_file_extension("test.ext1.bz"), "ext1")
        self.assertEqual(readers.get_file_extension("test.bz"), "")
        self.assertEqual(readers.get_file_extension(".bz"), "")
        self.assertEqual(readers.get_file_extension(""), "")

    def test_get_base_file_name(self):
        self.assertEqual(readers.get_base_file_name(test_files.ALIGNMENT_FILE_PATH_VALID_F),
                         "S1_11_S11_R1_001.fastqAligned.sortedByCoord.out")

        with self.assertRaises(IOError):
            readers.get_base_file_name("this.file.does.not.exist")
            readers.get_base_file_name(test_files.ALIGNMENT_FILE_PATH_D)


class BamReaderTest(unittest.TestCase):

    def test_check_file_validity_valid_f(self):
        with console_output() as (out, err):
            bamreader = BamReader(test_files.ALIGNMENT_FILE_PATH_VALID_F)
        self.assertEqual(bamreader.file_basename, "S1_11_S11_R1_001.fastqAligned.sortedByCoord.out")
        self.assertEqual(bamreader.compression, readers.COMPRESSION_None)
        self.assertEqual(bamreader.extension, bamreader.EXTENSION_BAM)
        self.assertIsNotNone(bamreader.alignment)
        self.assertIsNotNone(bamreader.alignment.bam_array)
        # FIXME cannot get stdout value after logging to streamhandler
        # self.assertTrue("Initialized BamReader" in out.getvalue())

    def test_check_file_validity_d(self):
        with self.assertRaises(IOError):
            check_file_path(test_files.ALIGNMENT_FILE_PATH_D)

    def test_check_file_validity_invalid_f(self):
        with self.assertRaises(IOError):
            check_file_path(test_files.ALIGNMENT_FILE_PATH_INVALID_F)

    def test_check_file_extension(self):
        with self.assertRaises(Exception):
            BamReader(test_files.ALIGNMENT_FILE_PATH_INVALID_EXT)


class FastaReaderTest(unittest.TestCase):

    def test_check_file_validity_valid_f(self):
        with console_output() as (out, err):
            fastareader = FastaReader(test_files.FASTA_FILE_PATH_VALID_F)
        self.assertEqual(fastareader.file_basename, "Bacillus_subtilis.ASM69118v1.dna.toplevel")
        self.assertEqual(fastareader.compression, readers.COMPRESSION_None)
        self.assertEqual(fastareader.extension, fastareader.EXTENSION_FA)
        # FIXME cannot get stdout value after logging to streamhandler
        # self.assertTrue("Initialized FastaReader" in out.getvalue())

    def test_check_file_validity_d(self):
        with self.assertRaises(IOError):
            check_file_path(test_files.FASTA_FILE_PATH_D)

    def test_check_file_validity_invalid_f(self):
        with self.assertRaises(IOError):
            check_file_path(test_files.FASTA_FILE_PATH_INVALID_F)

    def test_check_file_extension(self):
        with self.assertRaises(Exception):
            FastaReader(test_files.FASTA_FILE_PATH_INVALID_EXT)


class AnnotationReaderTest(unittest.TestCase):

    def test_check_file_validity_valid_f(self):
        with console_output() as (out, err):
            annotation_reader = AnnotationReader(test_files.ANNOT_FILE_PATH_VALID_F)
        self.assertEqual(annotation_reader.file_basename, "Bacillus_subtilis.ASM69118v1.37")
        self.assertEqual(annotation_reader.compression, readers.COMPRESSION_None)
        self.assertEqual(annotation_reader.extension, annotation_reader.EXTENSION_GFF3)

        self.assertIsNotNone(annotation_reader.annotation)
        self.assertIsNotNone(annotation_reader.annotation.transcript_assembly)
        # FIXME cannot get stdout value after logging to streamhandler
        # self.assertTrue("Initialized AnnotationReader" in out.getvalue())

    def test_check_file_validity_d(self):
        with self.assertRaises(IOError):
            check_file_path(test_files.FASTA_FILE_PATH_D)

    def test_check_file_validity_invalid_f(self):
        with self.assertRaises(IOError):
            check_file_path(test_files.FASTA_FILE_PATH_INVALID_F)

    def test_check_file_extension(self):
        with self.assertRaises(Exception):
            AnnotationReader(test_files.ANNOT_FILE_PATH_INVALID_EXT)


@contextmanager
def console_output():
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err
