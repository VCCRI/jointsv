import subprocess
import unittest


class TestJointSv(unittest.TestCase):
    log = True

    # Pre clean former results if any
    def setUp(self):
        subprocess.run(['rm', '-rf', 'testoutput'])

    # TODO We should refactor jointsv to allow us to output to stdout so that we don't need to create output test files
    # TODO We should provide more meaningful test failure output, for example a diff between the files
    def test_simple(self):
        self.assertEqual(1, 1, 'Should be 1')
        subprocess.run(['mkdir', 'testoutput'])
        processOutput = subprocess.run(
            ['python3', 'jointsv.py', 'test-data/simple/input/', '-o', 'testoutput/output.vcf'], capture_output=True)
        with open('testoutput/output.vcf') as expected_output_file:
            with open('test-data/simple/output/sample-output.vcf') as output_file:
                self.assertEqual(expected_output_file.readline(), output_file.readline(),
                                 "File Format (line 1) don't match")
                expected_headers = read_headers_from_file(expected_output_file)
                headers = read_headers_from_file(output_file)
                # Make sure all info and format headers are there regardless of order
                self.assertCountEqual(expected_headers, headers)
                expected_data_header = expected_headers[
                    len(expected_headers) - 1].split('\t')  # expected_output_file.readline().split("\t")
                data_header = headers[len(headers) - 1].split('\t')  # expected_output_file.readline().split("\t")
                # All data headers need to be there, samples may have swapped positions
                self.assertCountEqual(expected_data_header, data_header)
                number_of_samples = 5
                position_of_first_sample = 9
                sample_map = get_sample_mapping(number_of_samples, position_of_first_sample, expected_data_header,
                                                data_header)
                sample_name_map = get_sample_name_map(number_of_samples, position_of_first_sample, expected_data_header)

                data_line_1 = output_file.readline().split('\t')
                expected_data_line = expected_output_file.readline().split('\t')
                while len(data_line_1) > 0 and len(expected_data_line) > 0:
                    # Compare chromosome-position data
                    for i in range(position_of_first_sample):
                        logMessage('Compare ' + data_line_1[i] + " - " + expected_data_line[i])
                        self.assertEqual(data_line_1[i], expected_data_line[i])

                    # Compare sample specific data
                    for i in range(position_of_first_sample, position_of_first_sample + number_of_samples):
                        logMessage('Compare ' + data_line_1[sample_map[i]] + " - " + expected_data_line[i])
                        self.assertEqual(data_line_1[sample_map[i]],
                                         expected_data_line[i],
                                         "The values " + data_line_1[sample_map[i]] + "and " + expected_data_line[
                                             i] + "don't match for " + sample_name_map[i])
                    data_line_1 = output_file.readline().split('\t')
                    expected_data_line = expected_output_file.readline().split('\t')

    subprocess.run(['rm', '-rf', 'testoutput'])


# Get sample mappings. We take the output one as reference. In the final output, the order of the samples
# is not guaranteed, so we can't compare line by line. We'll use it as column mapping. If sample2 is in
# column 10 in expected output and in 14 in the output file, the mapping will be 10->14
def get_sample_mapping(number_of_samples, position_of_first_sample, expected_data_header, data_header):
    sample_map = {}
    current_sample = 0
    while current_sample < number_of_samples:
        logMessage("Finding..." + expected_data_header[position_of_first_sample + current_sample])
        position_of_sample_in_output_file = data_header.index(
            expected_data_header[position_of_first_sample + current_sample])
        sample_map[position_of_first_sample + current_sample] = position_of_sample_in_output_file
        current_sample += 1
    return sample_map


# Get sample names in the expected in a position-name map
def get_sample_name_map(number_of_samples, position_of_first_sample, expected_data_header):
    current_sample = 0
    # For debugging and logging purposes
    sample_names = {}
    while current_sample < number_of_samples:
        sample_names[position_of_first_sample + current_sample] = expected_data_header[
            position_of_first_sample + current_sample]
        current_sample += 1
    return sample_names


def logMessage(message):
    if TestJointSv.log:
        print(message)


def read_headers_from_file(file):
    header_list = []
    line = file.readline()
    still_in_headers = True
    while still_in_headers:
        header_list.append(line)
        if line.startswith("#CHROM"):
            still_in_headers = False
        else:
            line = file.readline()
    return header_list


if __name__ == '__main__':
    unittest.main()
