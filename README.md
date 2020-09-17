# JointSV

JointSV analyses multiple BND sample files in VCF format produced by Gridss
in order to detect Structural Variations by the Joint-call method.

# User manual

## Getting JointSV

Clone or download [the GitHub repository](https://github.com/VCCRI/jointsv).

## Requirements

JointSV requires Python 3.7 (newer versions may also work).

Some Python libraries are also required. They can be downloaded using the `pip` tool
from the directory that contains the JointSV source code.

```shell script
$ pip install -r requirements.txt
```

## Running the program

Place each input sample in a separate VCF file. JointSV does not support multi-sample
inputs.

Run the program like this:

```shell script
$ python3 jointsv.py input/*.vcf -o output.vcf
```

Alternatively, instead of enumerating the input files or using wildcards, a directory
can be used instead. In that case, all files (of any extension) in that directory will be
used as input. Subdirectories will not be traversed. This is useful to overcome the
limitations around the maximum number of files in the CLI. Example:

```shell script
$ python3 jointsv.py input/ -o output.vcf
```

## Command line options

* `FILE1, FILE2...`: name of the input files. If any of them is actually a directory,
  then all files in that directory will be used as input.
  The files will be processed in alphabetical filename order, regardless of the order
  in which they appear in the command line.

* `-o`/`--output FILE`: name of the output file. If omitted, the default is `output.vcf`.

* `--chromosome chr3,chr7`: indicates which chromosomes are to be processed, as comma-separated values.
  If this parameter is omitted, all chromosomes are included.

# Contact

Dr. Emma Rath - Victor Chang Cardiac Research Institute

### Collaborators

The following Atlassian employees have participated in writing this software
through a collaboration project between the Atlassian Foundation and the Victor
Chang Cardiac Research Institute:

* Luis Lainez
* Peter Grasevski
* Diego Berrueta
