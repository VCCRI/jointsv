import sys
from collections import defaultdict

from argument_parser import parse_arguments
from file_reader import read_records_from_files
from record_helper import *


def process_record_list(key, record_list, sample_names_to_header):
    # Create as many columns as samples
    # Process SVs if possible
    # if not possible return raw BNDs
    print("TODO: process record list at", key, len(record_list))
    creates_sv = False
    candidates = []
    for record in record_list:
        if is_trusted_record(record):
            creates_sv = compare_record_to_other_candidates(record, candidates)
        if creates_sv == True:
            break
        else:
            candidates.append(record)
    if (creates_sv == True):
        output = generate_sv_record(record_list)
    else:
        output = generate_non_sv_records(record_list, sample_names=sample_names_to_header.keys())
    return output


def compare_record_to_other_candidates(record, candidates):
    # TODO Some records will have the actual call in ALT, so we can't assume that mate_pos exists
    # We know that start position is the same and that the record is trusted
    for candidate_record in candidates:
        if are_pair_records(record, candidate_record):
            return True
    return False


def can_call_sv(key, record_list):
    higher_level_calls = ["DUP:INS", "DUP:TANDEM", "DUP:INS", "DEL", "INDEL"]
    # TODO
    can_call_sv = False
    for record in record_list:
        # If the record is called in its own merits in any sample <- Return true
        print(record)
        if record.ALT in higher_level_calls:
            return True

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

    return [] # TODO


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

    if original_record:
        call_data["GT"] = "0/1" # TODO: how to calculate this?
        call_data["TRANCHE2"] = original_record.INFO["TRANCHE2"]
        call_data["VAF"] = float(original_record.INFO["BNDVAF"])

    return vcfpy.Call(sample=sample_name, data=call_data)


def generate_non_sv_records(colocated_records, sample_names):
    """
    This function processes records that have not been used to call a SV.
    :param colocated_records:
    :param sample_names:
    :return:
    """

    # The colocated records need to be re-grouped based on their similary and actual position
    subgrouping_function = lambda record: (record.CHROM, record.POS, record.REF, str(record.ALT)) # TODO: excluded INFO because otherwise they don't match
    records_grouped_by_all_coordinates = group_by(colocated_records, key=subgrouping_function)

    # Once the regrouping has happened, each group will generate a single line in the output
    output = []
    for subkey, group in records_grouped_by_all_coordinates.items():
        print("Processing", subkey, group)

        # Build a map to easily find the records by the sample name
        sample_names_to_record = {record.ID[0]: record for record in group}

        # Generate calls for each sample in this group
        calls = [get_sample_call(sample_name, sample_names_to_record.get(sample_name, None))
                 for sample_name in sample_names]

        # Add a record to the output
        output.append(vcfpy.Record(
            CHROM=group[0].CHROM,  # by construction, all the grouped records have the same
            POS=group[0].POS,  # by construction, all the grouped records have the same
            ID=[group[0].CHROM + "_" + str(group[0].POS)],
            REF=group[0].REF,  # by construction, all the grouped records have the same
            ALT=group[0].ALT,  # by construction, all the grouped records have the same
            QUAL=None,  # FIXME: what to use here
            FILTER=[],  # FIXME: what to use here
            INFO=vcfpy.OrderedDict(),  # FIXME: what to use here
            FORMAT=["GT", "TRANCHE2", "VAF"],
            calls=calls))

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
    header = vcfpy.Header()
    header.add_format_line(vcfpy.OrderedDict(ID="GT", Number=1, Type="String", Description="Genotype"))
    header.add_format_line(vcfpy.OrderedDict(ID="TRANCHE2", Number=1, Type="String", Description="Quality category of GRIDSS structural variant calls determined using FILTER,SRQ,AS,RAS. Values are LOW INTERMEDIATE HIGH"))
    header.add_format_line(vcfpy.OrderedDict(ID="VAF", Number=1, Type="Float", Description="VAF of this SV call, derived from BNDVAF values of BND calls used to call this SV"))
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
