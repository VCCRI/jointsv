import sys
from argument_parser import parse_arguments
from file_reader import read_records_from_files
import vcfpy


def process_record_list(key, record_list):
    # Create as many columns as samples
    # Process SVs if possible
    # if not possible return raw BNDs
    print("TODO: process record list at", key, len(record_list))
    if can_call_sv(key, record_list):
        return generate_sv_record(record_list)
    else:
        return generate_non_sv_records(record_list)


def can_call_sv(key, record_list):
    # TODO
    return False


def generate_sv_record(record_list):
    # TODO
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


def generate_non_sv_records(record_list):
    # TODO
    return []


def write_output(output_list, output_file_path, sample_names):
    """
    Serialises the data into a VCF file.

    :param output_list:
    :param output_file_path:
    :return:
    """
    header = vcfpy.Header(lines=[],
                          samples=vcfpy.SamplesInfos(sample_names))

    writer = vcfpy.Writer.from_path(output_file_path, header)

    # Sort the output by chromosome and position
    sorting_function = lambda record: (record.CHROM, record.POS)
    output_list.sort(key=sorting_function)

    for output_record in output_list:
        writer.write_record(output_record)

    writer.close()


def main(args):
    print("Starting JointSV")
    args = parse_arguments(args)
    input_file_path_list = args.files
    output_file_path = args.output_file
    chromosome_set = set(args.chromosome_list) if args.chromosome_list else None

    # First, read all the data and group the records by CHROM + POS
    print("Reading inputs...")
    (records, sample_names) = read_records_from_files(input_file_path_list, chromosome_set)

    # Then process each group separately
    print("Processing", len(records), "groups from samples", sample_names, "...")
    output_list = []
    for key, colocated_records in records.items():
        output_list.extend(process_record_list(key, colocated_records))

    # Finally, write the output to a file
    print("Writing output...")
    write_output(output_list, output_file_path, sample_names)


if __name__ == "__main__":
    main(sys.argv)
