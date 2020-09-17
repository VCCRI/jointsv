import sys
from collections import defaultdict
from statistics import mean

from argument_parser import parse_arguments
from file_reader import read_records_from_files
from file_writer import write_output
from record_helper import *
from BndComparisonResult import BndComparisonResult
from sv_detector import *
import logging
import resource
import gc


def process_record_list(key, record_list, sample_names_to_header):
    # Create as many columns as samples
    # Process SVs if possible
    # if not possible return raw BNDs
    record_comparison_response = BndComparisonResult(False, None, None, None)
    candidates = []
    for record in record_list:
        if is_trusted_record(record):
            record_comparison_response = compare_record_to_other_candidates(record, candidates)
        if record_comparison_response.is_sv:
            break
        else:
            candidates.append(record)
    if record_comparison_response.is_sv:
        output = [
            generate_sv_record(record_list, record_comparison_response, sample_names=sample_names_to_header.keys())]
    else:
        output = generate_non_sv_records(record_list, sample_names=sample_names_to_header.keys())
    return output


def compare_record_to_other_candidates(record, candidates):
    # TODO Some records will have the actual call in ALT, so we can't assume that mate_pos exists
    # We know that start position is the same and that the record is trusted
    return_obj = BndComparisonResult(False, None, None, None)
    for candidate_record in candidates:
        if is_record_an_sv(record):
            return_obj.is_sv = True
            return_obj.type = get_alt_type(record)
            return_obj.initial_position = get_start_position(record)
            return_obj.final_position = get_end_position(record)
        if are_pair_records(record, candidate_record):
            return_obj.is_sv = True
            return_obj.type = extract_sv_type_from_record_pair(candidate_record, record)
            return_obj.initial_position = min(get_start_position(record), get_start_position(candidate_record))
            return_obj.final_position = max(get_end_position(record), get_end_position(candidate_record))
    return return_obj


def generate_sv_record(records, comparison_result, sample_names):
    """
    This method generates a single SV record after a call has been made over a set of input records
    :param records: the input records involved in the SV call
    :param comparison_result:
    :param sample_names:
    :return:
    """

    # Build a map to easily find the records by the sample name. It can be multi-valued
    sample_names_to_records = group_by(records, lambda record: get_sample_name(record))

    # Generate calls for each sample in this group
    calls = [get_sample_call(sample_name, sample_names_to_records.get(sample_name, None))
             for sample_name in sample_names]

    first_record_of_the_group = records[0]
    chrom = first_record_of_the_group.CHROM
    id_of_new_record = generate_id(chrom, comparison_result.initial_position)
    info = vcfpy.OrderedDict()
    info["SVTYPE"] = comparison_result.type
    info["END"] = comparison_result.final_position

    return vcfpy.Record(
        CHROM=chrom,  # by construction, all the grouped records have the same
        POS=comparison_result.initial_position,  # by construction, all the grouped records have the same
        ID=[id_of_new_record],
        REF=first_record_of_the_group.REF,  # by construction, all the grouped records have the same
        ALT=[vcfpy.Substitution(type_=comparison_result.type, value='<{}>'.format(comparison_result.type))],
        QUAL=None,  # FIXME: what to use here
        FILTER=[],  # FIXME: what to use here
        INFO=info,
        FORMAT=["GT", "TRANCHE2", "VAF"],
        calls=calls)


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


def maximum_tranche(records):
    tranches = set([get_tranche_2(record) for record in records])
    if "HIGH" in tranches:
        return "HIGH"
    elif "INTERMEDIATE" in tranches:
        return "INTERMEDIATE"
    elif "LOW" in tranches:
        return "LOW"
    else:
        return None


def get_sample_call(sample_name, records):
    """
    This function generates the Call for a single sample at at a given location, given a single record, multiple records or no record at all
    :param sample_name:
    :param records:
    :return:
    """
    call_data = vcfpy.OrderedDict.fromkeys(["GT", "TRANCHE2", "VAF"])

    if records:
        average_vaf = mean([float(record.INFO["BNDVAF"]) for record in records])
        call_data["GT"] = get_gt(average_vaf)
        call_data["TRANCHE2"] = maximum_tranche(records)
        call_data["VAF"] = average_vaf

    return vcfpy.Call(sample=sample_name, data=call_data)


def generate_id(chrom, pos):
    return chrom + "_" + str(pos)


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
        sample_names_to_record = group_by(group, get_sample_name)

        # Generate calls for each sample in this group
        calls = [get_sample_call(sample_name, sample_names_to_record.get(sample_name, []))
                 for sample_name in sample_names]

        # Add a record to the output
        first_record_of_the_group = group[0]
        id_of_new_record = generate_id(first_record_of_the_group.CHROM, first_record_of_the_group.POS)
        info = vcfpy.OrderedDict()
        if "END" in first_record_of_the_group.INFO:
            info["END"] = first_record_of_the_group.INFO[
                "END"]  # by construction, all the grouped records have the same
        if "INSSEQ" in first_record_of_the_group.INFO:
            info["INSSEQ"] = first_record_of_the_group.INFO[
                "INSSEQ"]  # by construction, all the grouped records have the same
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
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s')
    logging.info("Starting JointSV")
    args = parse_arguments(args)
    input_file_path_list = args.files
    output_file_path = args.output_file
    chromosome_set = set(args.chromosome_list) if args.chromosome_list else None

    log_resource_consumption()

    # First, read all the data and group the records by CHROM + POS
    (records, sample_names_to_header) = read_records_from_files(input_file_path_list, chromosome_set)
    sample_names = sample_names_to_header.keys()

    log_resource_consumption()

    # Then process each group separately
    logging.info("Processing %d co-located groups from %d samples", len(records), len(sample_names))
    output_records = []
    for key, colocated_records in records.items():
        output_records.extend(process_record_list(key, colocated_records, sample_names_to_header))

    log_resource_consumption()

    # Finally, write the output to a file
    logging.info("Writing %d records to output file '%s'", len(output_records), output_file_path)
    write_output(output_records, output_file_path, sample_names)

    log_resource_consumption()
    logging.info("JointSV finished successfully")


def log_resource_consumption():
    logging.debug("Resource usage: %s / GC generations: %s / GC stats: %s",
                  str(resource.getrusage(resource.RUSAGE_SELF)),
                  str(gc.get_count()),
                  str(gc.get_stats()))


if __name__ == "__main__":
    main(sys.argv)
