import re as regex
from record_helper import *

def extract_sv_type_from_record_pair(record1, record2):
    if is_del_sv(record1, record2):
        return "DEL"
    if is_a_dup_tandem_sv(record1, record2):
        return "DUP:TANDEM"
    return "UNKNOWN"


# DUP:TANDEM	CHROM	POS	ALT
# 		1	100	]1:200]N
#		1	200	N[1:100[
def is_a_dup_tandem_sv(record1, record2):

    first_record = record1 if get_start_position(record1) < get_start_position(record2) else record2
    second_record = record1 if get_start_position(record1) >= get_start_position(record2) else record2
    if (
            get_mate_orientation(first_record) == '-' and
            get_orientation(first_record) == '+' and
            len(get_sequence(first_record)) == 1 and
            get_mate_orientation(second_record) == '+' and
            get_orientation(second_record) == '-' and
            len(get_sequence(second_record)) == 1
    ):
        return True
    return False


# DEL: CHROM	POS	ALT
# 		1	100	N[1:200[
#		1	200	]1:100]N
def is_del_sv(record1, record2):
    first_record = record1 if get_start_position(record1) < get_start_position(record2) else record2
    second_record = record1 if get_start_position(record1) >= get_start_position(record2) else record2

    if (
            get_mate_orientation(first_record) == '+' and
            get_orientation(first_record) == '-' and
            len(get_sequence(first_record)) == 1 and
            get_mate_orientation(second_record) == '-' and
            get_orientation(second_record) == '+' and
            len(get_sequence(second_record)) == 1
    ):
        return True
    return False
