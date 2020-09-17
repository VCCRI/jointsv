import subprocess
import unittest


class TestJointSv(unittest.TestCase):

    # TODO We should refactor jointsv to allow us to output to stdout so that we don't need to create output test files
    # TODO We should provide more meaningful test failure output, for example a diff between the files
    def test_simple(self):
        self.assertEqual(1, 1, 'Should be 1')
        subprocess.run(['mkdir', 'testoutput'])
        subprocess.run(['python3', 'jointsv.py', 'test-data/simple/input/', '-o', 'testoutput/output.vcf'], capture_output=True)
        with open('testoutput/output.vcf') as f:
            output = f.read()
        with open('test-data/simple/output/sample-output.vcf') as f:
            expected_output = f.read()
        subprocess.run(['rm', '-rf', 'testoutput'])
        self.assertEqual(output, expected_output)

    # TODO We should add lower level tests that assert for a given input, we expect certain SVs to be called and
    #  exist in the output instead of just asserting against a whole output file

if __name__ == '__main__':
    unittest.main()
