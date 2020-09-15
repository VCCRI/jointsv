import sys
from argument_parser import *

from collections import defaultdict


def open_file(output_file_path):
    #open file
    print("TODO")

def read_records_from_files(filePath):
    # ...
    return defaultdict(list)

def process_record_list(record_list):
    # Create as many columns as samples
    # Process SVs if possible
    # if not possible return raw BNDs
    print("TODO")

def sort_elements(output_list):
    #sort the list before by start position
    return output_list

def write_element(item, file):
    #write element in file
    print("TODO")

def main(args):
    print("Starting Join SV")
    args = parse_arguments(args)
    filePath = args.files #read from stdin
    output_file_path = args.outputFile # get from stdin
    output_file = open_file(output_file_path)
    records = read_records_from_files(filePath)
    output_list = []
    for key,value in records.items():
        output_list.append(process_record_list(value))
    for item in sort_elements(output_file_path):
        write_element(item, output_file)


if __name__ == "__main__":
    main(sys.argv)
