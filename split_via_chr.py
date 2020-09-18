import os
def split_via_chr(list_of_files):
    chr_file_map = {}
    headers = ""
    for filepath in list_of_files:
        with open(filepath) as openfileobject:
            data_started = False
            for line in openfileobject:
                if not data_started:
                    headers += line
                    if line.startswith("#CHROM"):
                        data_started = True
                else:
                    chr = line.split('\t')[0]
                    if not chr in chr_file_map:
                        chr_file_map[chr] = open(chr + ".vcf", "w")
                        chr_file_map[chr].write(headers)
                    chr_file_map[chr].write(line)
    for chr in chr_file_map:
        chr_file_map[chr].flush()
        chr_file_map[chr].close()

def main():
    list_of_files = ["/Users/llainez/Documents/Development/jointsv/test-data/simple/input/sample1.vcf",
                     "/Users/llainez/Documents/Development/jointsv/test-data/simple/input/sample2.vcf",
                     "/Users/llainez/Documents/Development/jointsv/test-data/simple/input/sample3.vcf",
                     "/Users/llainez/Documents/Development/jointsv/test-data/simple/input/sample4.vcf",
                     "/Users/llainez/Documents/Development/jointsv/test-data/simple/input/sample5.vcf"]
    split_via_chr(list_of_files)

if __name__ == "__main__":
    main()
