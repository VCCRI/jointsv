from collections import defaultdict
from os import listdir
from os.path import isfile, join, basename
import vcfpy


def read_records_from_files(file_path_list, chromosome_set):
    multi_map = defaultdict(list)
    sample_names_to_header = {}
    for file_path in collect_all_file_names(file_path_list):
        reader = vcfpy.Reader.from_path(file_path)
        assert len(reader.header.samples.names) == 1, "Only records with exactly 1 sample are supported"
        sample_name = reader.header.samples.names[0]
        sample_names_to_header[sample_name] = reader.header
        for record in reader:
            if chromosome_set is None or record.CHROM in chromosome_set:
                assert len(record.ID) == 1, "Only records with exactly 1 ALT are supported"
                assert len(record.ALT) == 1, "Only records with exactly 1 ALT are supported"

                # Rewrite the ID to keep track of the sample where the data came from
                record.ID = (sample_name, record.ID[0])

                # Generate a key to group the records
                min_pos = min(record.POS, record.ALT[0].mate_pos) if isinstance(record.ALT, vcfpy.BreakEnd) else record.POS
                key = (record.CHROM, min_pos)

                # Put the record in a multimap grouped by its coordinates
                multi_map[key].append(record)
        reader.close()

    return multi_map, sample_names_to_header


def collect_all_file_names(file_path_list):
    files = []
    for file_path in file_path_list:
        if is_file_path_a_directory(file_path):
            files.extend(read_all_file_names_from_directory(file_path))
        else:
            files.append(file_path)
    return files


def is_file_path_a_directory(file_path):
    return file_path.endswith('/')


def read_all_file_names_from_directory(directory_path):
    return [directory_path + f for f in listdir(directory_path) if isfile(join(directory_path, f))]
