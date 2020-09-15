from collections import defaultdict
from os import listdir
from os.path import isfile, join, basename
import vcfpy


def read_records_from_files(file_path_list, chromosome_set):
    multi_map = defaultdict(list)
    sample_names = set()
    for file_path in collect_all_file_names(file_path_list):
        sample_name = basename(file_path)
        sample_names.add(sample_name)
        reader = vcfpy.Reader.from_path(file_path)
        for record in reader:
            if chromosome_set is None or record.CHROM in chromosome_set:
                assert len(record.ALT) == 1, "Only records with exactly 1 ALT are supported"

                # Rewrite the ID to keep track of the sample where the data came from
                record.ID = (sample_name, record.ID)

                # Generate a key to group the records
                key = (record.CHROM, min(record.POS, record.ALT[0].mate_pos))

                # Put the record in a multimap grouped by its coordinates
                multi_map[key].append(record)
        reader.close()

    return (multi_map, sample_names)


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
