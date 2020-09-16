import re as regex
from record_helper import *

re_left_N = regex.compile('\][a-zA-Z0-9.]+\:[0-9]+\][ACGTN]$')  # eg. ]1:200]N
re_left_NNN = regex.compile(
    '\][a-zA-Z0-9.]+\:[0-9]+\][ACGTN]{2,}$')  # eg. ]1:123200]TA, ]1:123200]TAC, ]1:123200]TACC, etc.

re_N_left = regex.compile('[ACGTN]\][a-zA-Z0-9.]+\:[0-9]+\]$')  # eg. T]1:123200]
re_NNN_left = regex.compile(
    '[ACGTN]{2,}\][a-zA-Z0-9.]+\:[0-9]+\]$')  # eg. TA]1:123200], TAC]1:123200], TACC]1:123200], etc.

re_right_N = regex.compile('\[[a-zA-Z0-9.]+\:[0-9]+\[[ACGTN]$')  # eg. [1:123200[T
re_right_NNN = regex.compile(
    '\[[a-zA-Z0-9.]+\:[0-9]+\[[ACGTN]{2,}$')  # eg. [1:123200[TA, [1:123200[TAC, [1:123200[TACC, etc.

re_N_right = regex.compile('[ACGTN]\[[a-zA-Z0-9.]+\:[0-9]+\[$')  # eg. N[1:100[
re_NNN_right = regex.compile(
    '[ACGTN]{2,}\[[a-zA-Z0-9.]+\:[0-9]+\[$')  # eg. TA[1:123200[, TAC[1:123200[, TACC[1:123200[, etc.


def extract_sv_type_from_record_pair(record1, record2):
    if is_a_dup_tandem_sv(record1, record2):
        return "ins"


# DUP:TANDEM	CHROM	POS	ALT
# 		1	100	]1:200]N
#		1	200	N[1:100[
def is_a_dup_tandem_sv(record1, record2):
    # 1 Check which one has the re_left_N and which one the re_N_right
    # 2 If the position of the one with re_left_N is less than that the other, is a DUP:TANDEM
    re_left_n_record = None
    re_n_right_record = None
    #if(record1.alt)
    return False
