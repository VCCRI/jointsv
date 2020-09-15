import sys
from argument_parser import parse_arguments
from file_reader import read_records_from_files
import vcfpy


def process_record_list(record_list):
    # Create as many columns as samples
    # Process SVs if possible
    # if not possible return raw BNDs
    print("TODO: process record list", len(record_list))
    format = vcfpy.OrderedDict.fromkeys(["sample1"], "format1")
    calls = [vcfpy.Call(sample="sample1", data=vcfpy.OrderedDict.fromkeys("key1", "value1")),
             vcfpy.Call(sample="sample2", data=vcfpy.OrderedDict.fromkeys("key2", "value2"))
             ]
    record = vcfpy.Record(CHROM="chr1",
                          POS=1000,
                          ID=["idFoo"],
                          REF="REF",
                          ALT=[vcfpy.SymbolicAllele("ALT")],
                          QUAL="QUAL",
                          FILTER=["FILTER"],
                          INFO=vcfpy.OrderedDict(),
                          FORMAT=format,
                          calls=calls)

    return [record]


def sort_elements(output_list):
    # sort the list before by start position
    return output_list


def write_output(output_list, output_file_path):
    """
    Serialises the data into a VCF file.

    :param output_list:
    :param output_file_path:
    :return:
    """
    header = vcfpy.Header(lines=[],
                          samples=vcfpy.SamplesInfos(["sample1", "sample2"]))

    writer = vcfpy.Writer.from_path(output_file_path, header)

    # Sort the output by chromosome and position
    sorting_function = lambda record: record.CHROM + "_" + str(record.POS)
    output_list.sort(key=sorting_function)

    for output_record in output_list:
        writer.write_record(output_record)

    writer.close()


def main(args):
    print("Starting JointSV")
    args = parse_arguments(args)
    file_path_list = args.files
    output_file_path = args.output_file
    chromosome_list = args.chromosome_list
    records = read_records_from_files(file_path_list, chromosome_list)
    output_list = []
    for key, value in records.items():
        output_list.extend(process_record_list(value))
    write_output(output_list, output_file_path)


if __name__ == "__main__":
    main(sys.argv)
