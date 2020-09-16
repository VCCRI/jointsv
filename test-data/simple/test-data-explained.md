## Simple test data

This data is designed to be the minimum input we could receive as input to the JointSV program that we can test the functionality with.

There are 5 pieces of sample input:
1. Has 2 BND records that can be called as a deletion SV on their own merit
2. Has 2 BND records that can be called as a deletion SV on their own merit, but with intermediate confidence
3. Has 2 BND records that are low confidence and therefore cannot be called on their own merit. However, they match with the records in #1 and #2 so should be called as SVs
4. Has 1 BND record and therefore cannot be called on its own merit, but since it matches with the existing SVs in sample #1 and #2 it should also be called as an SV
5. Has 2 BND records of low confidence that do not match any called SVs from other samples, therefore should not be called as an SV, and the original BND records should be outputted.

Based on this, the sample output should have:
* One called DEL SV, which matches samples 1, 2, 3 and 4. It does not match sample 5.
* Two remaining BND records from sample 5

### Simplification of data

In order to make these input files easy to read and understand, we've done the following:
* we've omitted values in the QUAL, FILTER, FORMAT, and SAMPLE columns, and replaced them with a dot (.). JointSV should not need to read the values in these fields for the input it receives, so it should function fine if that data is missing. In reality, these fields will contain content.
* we only include a small subset of values in the INFO field. In reality, there are dozens of values in this field.
* we only include information for a single chromosome (chr1). In reality, there can be hundreds of chromosome entries
