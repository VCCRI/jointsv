import sys
from argument_parser import parse_arguments
from file_reader import read_records_from_files
from collections import defaultdict
import vcfpy


def process_record_list(key, record_list, sample_names_to_header):
    # Create as many columns as samples
    # Process SVs if possible
    # if not possible return raw BNDs
    print("TODO: process record list at", key, len(record_list))
    if can_call_sv(key, record_list):
        return generate_sv_record(record_list)
    else:
        return generate_non_sv_records(record_list, sample_names_to_header)


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


def group_by(iterable, key):
    result = defaultdict(list)
    for item in iterable:
        result[key(item)].append(item)
    return result


def get_sample_call(sample_name, original_record):
    """
    This function generates the Call for a single sample at at a given location
    :param sample_name:
    :param original_record:
    :return:
    """
    call_data = vcfpy.OrderedDict.fromkeys(["GT", "TRANCHE2", "VAF"])
    # We may not have data in this group for all samples, so we fill the blanks with default values
    call_data.setdefault(".")
    if original_record:
        call_data["GT"] = "42" # TODO: how to calculate this?
        call_data["TRANCHE2"] = original_record.INFO["TRANCHE2"]
        call_data["VAF"] = original_record.INFO["BNDVAF"]

    return vcfpy.Call(sample=sample_name, data=call_data)


def generate_non_sv_records(record_list, sample_names_to_header):

    sample_names = sample_names_to_header.keys()
    subkey_func = lambda record: (record.CHROM, record.POS, record.REF, str(record.ALT)) # TODO: excluded INFO because otherwise they don't match
    format = ["GT", "TRANCHE2", "VAF"]
    output = []
    for subkey, group in group_by(record_list, key=subkey_func).items():
        print("Processing", subkey, group)
        sample_names_to_record = {record.ID[0]: record for record in group}

        calls = [get_sample_call(sample_name, sample_names_to_record.get(sample_name, None))
                 for sample_name in sample_names]

        record = vcfpy.Record(CHROM=group[0].CHROM, # by construction, all the grouped records have the same
                              POS=group[0].POS, # by construction, all the grouped records have the same
                              ID=[group[0].CHROM + "_" + str(group[0].POS)],
                              REF=group[0].REF, # by construction, all the grouped records have the same
                              ALT=group[0].ALT, # by construction, all the grouped records have the same
                              QUAL=None, # FIXME: what to use here
                              FILTER=[], # FIXME: what to use here
                              INFO=vcfpy.OrderedDict(), # FIXME: what to use here
                              FORMAT=format,
                              calls=calls)
        output.append(record)
    return output


def write_output(output_list, output_file_path, sample_names_to_header):
    """
    Serialises the data into a VCF file.

    :param output_list:
    :param output_file_path:
    :param sample_names_to_header:
    :return:
    """
    assert len(sample_names_to_header) > 0, "At least one sample is required"
    header = list(sample_names_to_header.values())[0].copy()
    header.samples = vcfpy.SamplesInfos(sample_names_to_header.keys())

    writer = vcfpy.Writer.from_path(output_file_path, header)

    # Sort the output by chromosome and position
    sorting_function = lambda record: (record.CHROM, record.POS)
    output_list.sort(key=sorting_function)

    for output_record in output_list:
        # Convert the IDs back to strings, because tuples cannot be serialised by VCFPy
        # output_record.ID = [str(output_record.ID[0]) + "_" + str(output_record.ID[1])]
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
    (records, sample_names_to_header) = read_records_from_files(input_file_path_list, chromosome_set)

    # Then process each group separately
    print("Processing", len(records), "groups from samples", sample_names_to_header.keys(), "...")
    output_list = []
    for key, colocated_records in records.items():
        output_list.extend(process_record_list(key, colocated_records, sample_names_to_header))

    # Finally, write the output to a file
    print("Writing output...")
    write_output(output_list, output_file_path, sample_names_to_header)


if __name__ == "__main__":
    main(sys.argv)
