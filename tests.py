import subprocess
import unittest


class TestJointSv(unittest.TestCase):

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


if __name__ == '__main__':
    unittest.main()
