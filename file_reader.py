from collections import defaultdict
from os import listdir
from os.path import isdir, isfile, join
from record_helper import *
import vcfpy
import logging
import gc
import io


class ChromosomeFilter(io.TextIOWrapper):
    """
    This class wraps a stream and modifies the behaviour of readline to only return the lines that are headers or
    are related to one of the chromosomes in the given set. This is done to prevent VCFPy from doing a heavy parsing
    of all the fields for the records that belong to irrelevant chromosomes.
    """
    def __init__(self, buffer, chromosome_set):
        io.TextIOWrapper.__init__(self, buffer)
        self.chromosome_set = chromosome_set

    def readline(self):
        line = self.buffer.readline()

        # Return immediately if no filtering is to be applied
        if not self.chromosome_set:
            return line

        while line:
            # Headers are allowed
            if line[0] == '#':
                return line

            # Check only the first column
            fields = line.split('\t')
            if len(fields) >= 8 and fields[0] in self.chromosome_set:
                return line

            # Continue processing the file
            line = self.buffer.readline()
            
        return line


def read_records_from_files(file_path_list, chromosome_set=None):
    """
    This method parses the provided input files, filters the selected chromosomes, and returns the records grouped by
    their location (where the location is not necessarily CHROM+POS, but it could be CHROM+END if the END is lower
    than the POS).
    The ID of each record if modified to make it a tuple (sample,original_id), in order to avoid clashes and also
    be able to trace each record to their sample.
    :param file_path_list: file paths to be read. If any of them is a directory, the files in the directory will be parsed (non recursively)
    :param chromosome_set: the set of chromosomes of interest. Use None to read all of them
    :return: a multimap where the key is the location (see note above), and the value is a list of all the records found at that location
    """
    colocated_records_multimap = defaultdict(list)
    all_file_paths = collect_all_file_names(file_path_list)
    logging.info("Reading %d input files", len(all_file_paths))

    record_count = 0
    sample_names = []
    for file_path in all_file_paths:
        gc.collect()
        logging.debug("Reading file '%s'", file_path)
        with open(file_path, "rt") as file:
            with vcfpy.Reader.from_stream(ChromosomeFilter(file, chromosome_set)) as reader:

                # Extract the sample name from the headers
                assert len(reader.header.samples.names) == 1, "Only records with exactly 1 sample are supported"
                sample_name = reader.header.samples.names[0]
                sample_names.append(sample_name)

                # Process each record in the file. Note that the VCFPy parser is a streaming parser. It does not hold
                # the whole AST in memory, instead it goes record by record.
                for record in reader:
                    # Improvement idea: instead of filtering by chromosome after the whole record has been parsed,
                    # the filtering could happen before parsing the whole record. Of course that can be accomplished
                    # with UNIX pipes, but it would be great to have a full Python solution.
                    if chromosome_set is None or record.CHROM in chromosome_set:
                        assert len(record.ID) == 1, "Only records with exactly 1 ID are supported"
                        assert len(record.ALT) == 1, "Only records with exactly 1 ALT are supported"

                        # Rewrite the ID to keep track of the sample where the data came from
                        record.ID = (sample_name, record.ID[0])

                        # Identify the location of the record so they can be grouped by co-located records. Note the
                        # position is not necessarily the one in the input BND, but it could be the end position if
                        # it references a previous position.
                        min_pos = get_start_position(record)
                        end_position = get_end_position(record)
                        if end_position is not None:
                            min_pos = min(get_start_position(record), get_end_position(record))
                        record_location = (record.CHROM, min_pos)

                        # Put the record in a multimap grouped by its coordinates
                        colocated_records_multimap[record_location].append(record)
                        record_count += 1

    logging.debug("Found %d records from %d samples and grouped them in %d co-located groups",
                  record_count, len(sample_names), len(colocated_records_multimap))

    return colocated_records_multimap, sample_names


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
