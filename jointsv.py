import sys
from collections import defaultdict

from argument_parser import parse_arguments
from file_reader import read_records_from_files
from file_writer import write_output
from record_helper import *
import logging


def process_record_list(key, record_list, sample_names_to_header):
    # Create as many columns as samples
    # Process SVs if possible
    # if not possible return raw BNDs
    logging.info("TODO: process %d records list at %s", len(record_list), str(key))
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


def get_gt(original_bndvat):
    if original_bndvat > 0.85:
        return "1/1"
    elif original_bndvat < 0.15:
        return "0/0"
    else:
        return "0/1"


def get_sample_call(sample_name, original_record):
    """
    This function generates the Call for a single sample at at a given location
    :param sample_name:
    :param original_record:
    :return:
    """
    call_data = vcfpy.OrderedDict.fromkeys(["GT", "TRANCHE2", "VAF"])

    if original_record:
        original_bndvaf = float(original_record.INFO["BNDVAF"])
        call_data["GT"] = get_gt(original_bndvaf)
        call_data["TRANCHE2"] = get_tranche_2(original_record)
        call_data["VAF"] = original_bndvaf

    return vcfpy.Call(sample=sample_name, data=call_data)


def generate_non_sv_records(colocated_records, sample_names):
    """
    This function processes records that have not been used to call a SV.
    :param colocated_records:
    :param sample_names:
    :return:
    """

    # The co-located records need to be re-grouped based not just on their true position (CHROM+POS) but also similarity
    subgrouping_function = lambda record: (record.CHROM,
                                           record.POS,
                                           record.REF,
                                           str(record.ALT),
                                           record.INFO.get("END", None),
                                           record.INFO.get("INSSEQ", None))
    records_grouped_by_all_coordinates = group_by(colocated_records, key=subgrouping_function)

    # Once the regrouping has happened, each group will generate exactly one line in the output. These lines
    # may be produced out-of-order, but we don't care because we will sort them later before generating the VCF.
    output = []
    for subkey, group in records_grouped_by_all_coordinates.items():
        # Build a map to easily find the records by the sample name
        sample_names_to_record = {get_sample_name(record): record for record in group}

        # Generate calls for each sample in this group
        calls = [get_sample_call(sample_name, sample_names_to_record.get(sample_name, None))
                 for sample_name in sample_names]

        # Add a record to the output
        first_record_of_the_group = group[0]
        id_of_new_record = first_record_of_the_group.CHROM + "_" + str(first_record_of_the_group.POS)
        info = vcfpy.OrderedDict()
        if "END" in first_record_of_the_group.INFO:
            info["END"] = first_record_of_the_group.INFO["END"] # by construction, all the grouped records have the same
        if "INSSEQ" in first_record_of_the_group.INFO:
            info["INSSEQ"] = first_record_of_the_group.INFO["INSSEQ"] # by construction, all the grouped records have the same
        output.append(vcfpy.Record(
            CHROM=first_record_of_the_group.CHROM,  # by construction, all the grouped records have the same
            POS=first_record_of_the_group.POS,  # by construction, all the grouped records have the same
            ID=[id_of_new_record],
            REF=first_record_of_the_group.REF,  # by construction, all the grouped records have the same
            ALT=first_record_of_the_group.ALT,  # by construction, all the grouped records have the same
            QUAL=None,  # FIXME: what to use here
            FILTER=[],  # FIXME: what to use here
            INFO=info,
            FORMAT=["GT", "TRANCHE2", "VAF"],
            calls=calls))

    return output


def main(args):
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')
    logging.info("Starting JointSV")
    args = parse_arguments(args)
    input_file_path_list = args.files
    output_file_path = args.output_file
    chromosome_set = set(args.chromosome_list) if args.chromosome_list else None

    # First, read all the data and group the records by CHROM + POS
    logging.info("Reading %d input files", len(input_file_path_list))
    (records, sample_names_to_header) = read_records_from_files(input_file_path_list, chromosome_set)
    logging.info("%d samples loaded", len(sample_names_to_header.keys()))

    # Then process each group separately
    logging.info("Processing %d groups from samples %s", len(records), sample_names_to_header.keys())
    output_list = []
    for key, colocated_records in records.items():
        output_list.extend(process_record_list(key, colocated_records, sample_names_to_header))

    # Finally, write the output to a file
    logging.info("Writing output file at %s", output_file_path)
    write_output(output_list, output_file_path, sample_names_to_header)
    logging.info("JointSV finished successfully")


if __name__ == "__main__":
    main(sys.argv)
