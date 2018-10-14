import unittest

from fivepseq.util.formatting import pad_spaces


class StringFormattingTest(unittest.TestCase):

    def test_pad_spaces(self):
        test_string1 = "a" * 10
        test_string2 = "b" * 20
        test_string3 = "c" * 30

        self.assertEqual(len(pad_spaces(test_string1)), 20)
        self.assertEqual(len(pad_spaces(test_string1, 30)), 30)
        self.assertEqual(len(pad_spaces(test_string2)), 20)
        self.assertEqual(len(pad_spaces(test_string3)), 30)
        self.assertEqual(len(pad_spaces(test_string1, -6)), 10)
