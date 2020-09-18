import sys

from argument_parser import parse_arguments
from file_reader import read_records_from_files
from file_writer import write_output
from record_helper import *
from sv_detector import *
from generate import generate_sv_record, generate_non_sv_records
from BndComparisonResult import BndComparisonResult
import logging
import resource
import gc


def process_record_list(key, record_list, sample_names):
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
        output = [generate_sv_record(record_list, record_comparison_response, sample_names)]
    else:
        output = generate_non_sv_records(record_list, sample_names)

    return output, record_comparison_response


def compare_record_to_other_candidates(record, candidates):
    # We know that start position is the same and that the record is trusted
    return_obj = BndComparisonResult(False, None, None, None)
    if is_record_an_sv(record):
        return_obj.is_sv = True
        return_obj.svtype = get_sv_type_from_record(record)
        return_obj.initial_position = get_start_position(record)
        return_obj.final_position = get_end_position(record)
        return_obj.insseq = get_insseq_from_sv(record)
    else:
        for candidate_record in candidates:
            if are_pair_records(record, candidate_record):
                return_obj.is_sv = True
                return_obj.svtype = extract_sv_type_from_record_pair(candidate_record, record)
                return_obj.initial_position = min(get_start_position(record), get_start_position(candidate_record))
                return_obj.final_position = max(get_end_position(record), get_end_position(candidate_record))
                return_obj.insseq = get_insseq_from_bnds(return_obj.svtype, candidate_record, record)
    return return_obj


def main(args):
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s')
    args = parse_arguments(args)
    input_file_path_list = args.files
    output_file_path = args.output_file
    chromosome_set = set(args.chromosome_list) if args.chromosome_list else None

    logging.info("Starting JointSV")
    jointsv(input_file_path_list, output_file_path, chromosome_set)
    logging.info("JointSV finished successfully")


def jointsv(input_file_path_list, output_file_path, chromosome_set):
    log_resource_consumption()

    # First, read all the data and group the records by CHROM + POS
    (records, sample_name_to_header) = read_records_from_files(input_file_path_list, chromosome_set)
    sample_names = sample_name_to_header.keys()

    log_resource_consumption()

    # Then process each group separately
    logging.info("Processing %d co-located groups from %d samples", len(records), len(sample_names))
    positions = list(records.keys())
    record_accumulator = []
    sv_calls = 0
    for position in positions:
        # Instead of iterating over the .items() of the dictionary, we make a copy of the list of keys and iterate
        # over it. This makes it possible to delete entries from the dictionary as we go, reducing the memory
        # footprint.
        colocated_records = records[position]
        output_records, record_comparison_response = process_record_list(position, colocated_records, sample_names)
        record_accumulator.extend(output_records)
        del records[position]
        if record_comparison_response.is_sv:
            sv_calls += 1

    log_resource_consumption()

    # Finally, write the output to a file
    logging.info("Writing %d records to output file '%s', including %d SV calls",
                 len(record_accumulator),
                 output_file_path,
                 sv_calls)
    write_output(record_accumulator, output_file_path, sample_name_to_header, chromosome_set)

    log_resource_consumption()


def log_resource_consumption():
    logging.debug("Resource usage: %s / GC generations: %s / GC stats: %s",
                  str(resource.getrusage(resource.RUSAGE_SELF)),
                  str(gc.get_count()),
                  str(gc.get_stats()))


if __name__ == "__main__":
    main(sys.argv)
