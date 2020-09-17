import subprocess
import unittest


class TestJointSv(unittest.TestCase):

    #Pre clean former results if any
    def setUp(self):
        subprocess.run(['rm', '-rf', 'testoutput'])
    # TODO We should refactor jointsv to allow us to output to stdout so that we don't need to create output test files
    # TODO We should provide more meaningful test failure output, for example a diff between the files
    def test_simple(self):
        self.assertEqual(1, 1, 'Should be 1')
        subprocess.run(['mkdir', 'testoutput'])
        processOutput = subprocess.run(['python3', 'jointsv.py', 'test-data/simple/input/', '-o', 'testoutput/output.vcf'], capture_output=True)
        with open('testoutput/output.vcf') as expected_output_file:
            with open('test-data/simple/output/sample-output.vcf') as output_file:
                self.assertEqual(expected_output_file.readline(),output_file.readline(), "File Format (line 1) don't match")
                expected_headers = read_headers_from_file(expected_output_file)
                headers = read_headers_from_file(output_file)
                self.assertCountEqual(expected_headers, headers) #assertListEquals regardless order


    subprocess.run(['rm', '-rf', 'testoutput'])

    # TODO We should add lower level tests that assert for a given input, we expect certain SVs to be called and
    #  exist in the output instead of just asserting against a whole output file

def read_headers_from_file(file):
    header_list = []
    line = file.readline()
    while(line.startswith("#")):
        header_list.append(line)
        line = file.readline()
    return header_list

if __name__ == '__main__':
    unittest.main()
