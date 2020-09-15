from collections import defaultdict


def read_records_from_files(file_path_list, chromosome_list):
    multi_map = defaultdict(list)
    for file_path in file_path_list:
        data_started = False
        file = open(file_path, 'r')
        while not data_started:
            line = file.readline()
            if line.startswith("#CHROM"):
                data_started = True
                

    return multi_map
