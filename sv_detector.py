from record_helper import *


def is_record_an_sv(record):
    return get_alt_type(record) != "BND"


def are_pair_records(record1, record2):
    return (
            get_end_position(record1) == get_start_position(record2) and
            get_start_position(record1) == get_end_position(record2) and
            get_sample_name(record1) == get_sample_name(record2)
    )


def extract_sv_type_from_record_pair(record1, record2):
    if is_del_sv(record1, record2):
        return "DEL"
    if is_a_dup_tandem_sv(record1, record2):
        return "DUP:TANDEM"
    if is_ins_sv(record1, record2):
        return "INS"
    if is_indel_sv(record1, record2):
        return "INDEL"
    if is_dup_ins_sv(record1, record2):
        return "DUP:INS"
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


# INS:		CHROM	POS	ALT
# 		1	100	]1:100]NNN
#		1	100	NNN[1:100[
def is_ins_sv(record1, record2):
    first_record = record1
    second_record = record2
    # They have the same POS and mate_pos
    if get_start_position(record1) != get_start_position(record2):
        return False
    if (
            get_mate_orientation(first_record) == '+' and
            get_orientation(first_record) == '-' and
            len(get_sequence(first_record)) > 1 and
            get_mate_orientation(second_record) == '-' and
            get_orientation(second_record) == '+' and
            len(get_sequence(second_record)) > 1
    ):
        return True
    elif (
            get_mate_orientation(second_record) == '+' and
            get_orientation(second_record) == '-' and
            len(get_sequence(second_record)) > 1 and
            get_mate_orientation(first_record) == '-' and
            get_orientation(first_record) == '+' and
            len(get_sequence(first_record)) > 1
    ):
        return True
    return False


# INDEL:	CHROM	POS	ALT
# 		1	100	NNN[1:200[
#		1	200	]1:100]NNN
def is_indel_sv(record1, record2):
    first_record = record1 if get_start_position(record1) < get_start_position(record2) else record2
    second_record = record1 if get_start_position(record1) >= get_start_position(record2) else record2

    if (
            get_mate_orientation(first_record) == '+' and
            get_orientation(first_record) == '-' and
            len(get_sequence(first_record)) > 1 and
            get_mate_orientation(second_record) == '-' and
            get_orientation(second_record) == '+' and
            len(get_sequence(second_record)) > 1
    ):
        return True
    return False


# DUP:INS:	CHROM	POS	ALT
# 		1	100	]1:200]NNN
#		1	200	NNN[1:100[
def is_dup_ins_sv(record1, record2):
    first_record = record1 if get_start_position(record1) < get_start_position(record2) else record2
    second_record = record1 if get_start_position(record1) >= get_start_position(record2) else record2

    if (
            get_mate_orientation(first_record) == '-' and
            get_orientation(first_record) == '+' and
            len(get_sequence(first_record)) > 1 and
            get_mate_orientation(second_record) == '+' and
            get_orientation(second_record) == '-' and
            len(get_sequence(second_record)) > 1
    ):
        return True
    return False


def get_insseq_from_sv(record):
    if "INSSEQ" in record.INFO:
        return record.INFO["INSSEQ"]
    return None

'''
The INSSEQ is the last 2 characters in the record with ALT with this shape:
    ]1:200]NNN
We will also append the REF if the position in ALT is the same as the position of the record
    #CHROM	POS	ALT REF
    #1  100 ABC[1:200[ D -> INSSEQ = BC
    #CHROM	POS	ALT REF
    #2  100 ABC[1:100[ D -> INSSEQ = BCD
'''
def get_insseq_from_bnds(sv_type, record1, record2):
    if sv_type in get_sv_type_that_require_insseq():
        record_containing_the_insseq = record1
        if record2.ALT[0].orientation == '-':
            record_containing_the_insseq = record2
        insseq = get_sequence(record_containing_the_insseq)
        # Very defensive programming
        if insseq is not None and len(insseq) > 1:
            insseq = insseq[1:]
            if get_start_position(record_containing_the_insseq) == get_end_position(record_containing_the_insseq):
                insseq += get_ref(record_containing_the_insseq)
        return insseq
    return None


def get_sv_type_that_require_insseq():
    return ["INS", "INDEL", "DUP:INS"]
