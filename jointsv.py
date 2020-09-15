import sys
from argument_parser import parse_arguments
from file_reader import read_records_from_files


def open_file(output_file_path):
    # open file
    print("TODO")


def process_record_list(record_list):
    # Create as many columns as samples
    # Process SVs if possible
    # if not possible return raw BNDs
    print("TODO")


def sort_elements(output_list):
    # sort the list before by start position
    return output_list


def write_element(item, file):
    # write element in file
    print("TODO")


def main(args):
    print("Starting Join SV")
    args = parse_arguments(args)
    file_path_list = args.files
    output_file_path = args.output_file
    chromosome_list = args.chromosome_list
    records = read_records_from_files(file_path_list, chromosome_list)
    output_file = open_file(output_file_path)
    output_list = []
    for key, value in records.items():
        output_list.append(process_record_list(value))
    for item in sort_elements(output_file_path):
        write_element(item, output_file)


if __name__ == "__main__":
    main(sys.argv)
