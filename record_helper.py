"""
This is a collection of methods that provides a friendly and streamlined interface to operate with the Record object.
Please refrain from operating with the record object itself, it's best to create a method that does that and use the
method itself, since we could change the behavior of the methods independently and these changes would be 'broadcasted'
throughout the software.
"""
def get_tranche_2(record):
    return record.INFO["TRANCHE2"]


def get_sample_name(record):
    return record.ID[0]


def get_start_position(record):
    return record.POS


def get_end_position(record):
    if hasattr(record.ALT[0], 'mate_pos'):
        return record.ALT[0].mate_pos
    if "END" in record.INFO:
        return record.INFO["END"]
    return None


def is_trusted_record(record):
    return get_tranche_2(record) == "HIGH" or get_tranche_2(record) == "INTERMEDIATE"


def get_orientation(record):
    return record.ALT[0].orientation


def get_mate_orientation(record):
    return record.ALT[0].mate_orientation


def get_sequence(record):
    return record.ALT[0].sequence


def get_alt_type(record):
    return record.ALT[0].type


def get_sv_type_from_record(record):
    if 'SVTYPE' in record.INFO:
        return record.INFO['SVTYPE']
    if hasattr(record.ALT[0], 'value'):
        return record.ALT[0].value
    return type


def get_ref(record):
    return record.REF
