from collections import defaultdict
from os import listdir
from os.path import isdir, isfile, join
from record_helper import *
import vcfpy
import logging
import gc


def read_records_from_files(file_path_list, chromosome_set):
    colocated_records_multimap = defaultdict(list)
    all_file_paths = collect_all_file_names(file_path_list)
    logging.info("Reading %d input files", len(all_file_paths))

    record_count = 0
    sample_names_to_header = {}
    for file_path in all_file_paths:
        gc.collect()
        logging.debug("Reading file '%s'", file_path)
        with vcfpy.Reader.from_path(file_path) as reader:

            # Extract the sample name from the headers
            assert len(reader.header.samples.names) == 1, "Only records with exactly 1 sample are supported"
            sample_name = reader.header.samples.names[0]
            sample_names_to_header[sample_name] = reader.header

            # Process each record in the file
            for record in reader:
                if chromosome_set is None or record.CHROM in chromosome_set:
                    assert len(record.ID) == 1, "Only records with exactly 1 ID are supported"
                    assert len(record.ALT) == 1, "Only records with exactly 1 ALT are supported"

                    # Rewrite the ID to keep track of the sample where the data came from
                    record.ID = (sample_name, record.ID[0])

                    # Identify the location of the record so they can be grouped by co-located records. Note the
                    # position is not necessarily the one in the input BND, but it could be the end position if
                    # it references a previous position.
                    min_pos = get_start_position(record)
                    if get_end_position(record):
                        min_pos = min(get_start_position(record), get_end_position(record))
                    record_location = (record.CHROM, min_pos)

                    # Put the record in a multimap grouped by its coordinates
                    colocated_records_multimap[record_location].append(record)
                    record_count += 1

    logging.debug("Found %d records from %d samples and grouped them in %d co-located groups",
                  record_count, len(sample_names_to_header), len(colocated_records_multimap))

    return colocated_records_multimap, sample_names_to_header


def collect_all_file_names(file_path_list):
    files = []
    for file_path in file_path_list:
        if isdir(file_path):
            files.extend(read_all_file_names_from_directory(file_path))
        else:
            files.append(file_path)
    return files


def read_all_file_names_from_directory(directory_path):
    return [join(directory_path, f) for f in listdir(directory_path) if isfile(join(directory_path, f))]
