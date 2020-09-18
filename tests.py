import unittest
import io
from contextlib import redirect_stdout
from jointsv import jointsv


class TestJointSv(unittest.TestCase):
    def test_simple(self):
        with io.StringIO() as actual_output_file:
            with redirect_stdout(actual_output_file):
                jointsv(['test-data/simple/input/'], output_file_path='-', chromosome_set=None)
            actual_lines = actual_output_file.getvalue().split('\n')

        with open('test-data/simple/output/sample-output.vcf') as expected_output_file:
            expected_lines = [line.rstrip('\n') for line in expected_output_file.readlines()]

        self.assertEqual(actual_lines[0], expected_lines[0], "File Format (line 1) don't match")

        actual_headers = filter_headers(actual_lines)
        expected_headers = filter_headers(expected_lines)

        # Make sure all info and format headers are there regardless of order
        self.assertCountEqual(actual_headers, expected_headers)

        # All data headers need to be there
        actual_data_header = actual_headers[-1].split('\t')
        expected_data_header = expected_headers[-1].split('\t')
        self.assertCountEqual(actual_data_header, expected_data_header)

        expected_data_lines = filter_data_lines(expected_lines)
        actual_data_lines = filter_data_lines(actual_lines)
        self.assertEqual(len(expected_data_lines), len(actual_data_lines))

        for expected_data_line, actual_data_line in zip(expected_data_lines, actual_data_lines):
            expected_fields = expected_data_line.split('\t')
            actual_fields = actual_data_line.split('\t')
            self.assertListEqual(expected_fields, actual_fields)


def filter_headers(lines):
    return [line for line in lines if line.startswith("##")]


def filter_data_lines(lines):
    return [line for line in lines if not line.startswith("##") and line != '']


if __name__ == '__main__':
    unittest.main()
