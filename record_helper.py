import vcfpy

def get_tranche_2(record):
    return record.INFO["TRANCHE2"]


def get_sample_name(record):
    return record.ID[0]


def get_start_position(record):
    return record.POS


def get_end_position(record):
    return record.ALT[0].mate_pos


def is_trusted_record(record):
    return get_tranche_2(record) == "HIGH" or get_tranche_2(record) == "INTERMEDIATE"


# Maybe name doesnt match the logic, but it's what we found so far
def is_record_an_sv(record):
    return isinstance(record.ALT, vcfpy.BreakEnd)


def are_pair_records(record1, record2):
    return (
            get_end_position(record1) == get_start_position(record2) and
            get_start_position(record1) == get_end_position(record2) and
            get_sample_name(record1) == get_sample_name(record2)
    )

def get_orientation(record):
    return record.ALT[0].orientation

def get_mate_orientation(record):
    return record.ALT[0].mate_orientation

def get_sequence(record):
    return record.ALT[0].sequence