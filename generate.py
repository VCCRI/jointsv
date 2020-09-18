from collections import defaultdict
from statistics import mean
from record_helper import *
import vcfpy


def generate_sv_record(records, comparison_result, sample_names):
    """
    This method generates a single SV record after a call has been made over a set of input records
    :param records: the input records involved in the SV call
    :param comparison_result:
    :param sample_names:
    :return:
    """

    # Build a map to easily find the records by the sample name. It can be multi-valued
    sample_names_to_records = group_by(records, lambda record: get_sample_name(record))

    # Generate calls for each sample in this group
    calls = [get_sample_call(sample_name, sample_names_to_records.get(sample_name, None))
             for sample_name in sample_names]

    first_record_of_the_group = records[0]
    chrom = first_record_of_the_group.CHROM
    id_of_new_record = generate_id(chrom, comparison_result.initial_position)
    info = vcfpy.OrderedDict()
    info["SVTYPE"] = comparison_result.svtype
    info["END"] = comparison_result.final_position
    if comparison_result.insseq is not None:
        info["INSSEQ"] = comparison_result.insseq
    return vcfpy.Record(
        CHROM=chrom,  # by construction, all the grouped records have the same
        POS=comparison_result.initial_position,  # by construction, all the grouped records have the same
        ID=[id_of_new_record],
        REF=first_record_of_the_group.REF,  # by construction, all the grouped records have the same
        ALT=[vcfpy.Substitution(type_=comparison_result.svtype, value='<{}>'.format(comparison_result.svtype))],
        QUAL=maximum_qual(records),
        FILTER=["PASS"],
        INFO=info,
        FORMAT=["GT", "TRANCHE2", "VAF"],
        calls=calls)


def generate_non_sv_records(colocated_records, sample_names):
    """
    This function processes records that have not been used to call a SV.
    :param colocated_records:
    :param sample_names:
    :return:
    """

    # The co-located records need to be re-grouped based not just on their true position (CHROM+POS) but also similarity
    subgrouping_function = lambda record: (record.CHROM,
                                           record.POS,
                                           record.REF,
                                           str(record.ALT),
                                           record.INFO.get("END", None),
                                           record.INFO.get("INSSEQ", None))
    records_grouped_by_all_coordinates = group_by(colocated_records, key=subgrouping_function)

    # Once the regrouping has happened, each group will generate exactly one line in the output. These lines
    # may be produced out-of-order, but we don't care because we will sort them later before generating the VCF.
    output = []
    for subkey, group in records_grouped_by_all_coordinates.items():
        # Build a map to easily find the records by the sample name
        sample_names_to_record = group_by(group, get_sample_name)

        # Generate calls for each sample in this group
        calls = [get_sample_call(sample_name, sample_names_to_record.get(sample_name, []))
                 for sample_name in sample_names]

        # Add a record to the output
        first_record_of_the_group = group[0]
        id_of_new_record = generate_id(first_record_of_the_group.CHROM, first_record_of_the_group.POS)
        info = vcfpy.OrderedDict()
        info["SVTYPE"] = "BND"
        info["TRANCHE2"] = maximum_tranche(group)
        info["BNDVAF"] = get_average_vaf(group)
        if "END" in first_record_of_the_group.INFO:
            # by construction, all the grouped records have the same
            info["END"] = first_record_of_the_group.INFO["END"]
        if "INSSEQ" in first_record_of_the_group.INFO:
            # by construction, all the grouped records have the same
            info["INSSEQ"] = first_record_of_the_group.INFO["INSSEQ"]
        output.append(vcfpy.Record(
            CHROM=first_record_of_the_group.CHROM,  # by construction, all the grouped records have the same
            POS=first_record_of_the_group.POS,  # by construction, all the grouped records have the same
            ID=[id_of_new_record],
            REF=first_record_of_the_group.REF,  # by construction, all the grouped records have the same
            ALT=first_record_of_the_group.ALT,  # by construction, all the grouped records have the same
            QUAL=maximum_qual(group),
            FILTER=["PASS"],
            INFO=info,
            FORMAT=["GT", "TRANCHE2", "VAF"],
            calls=calls))

    return output


def group_by(iterable, key):
    result = defaultdict(list)
    for item in iterable:
        result[key(item)].append(item)
    return result


def get_gt(original_bndvat):
    if original_bndvat > 0.85:
        return "1/1"
    elif original_bndvat < 0.15:
        return "0/0"
    else:
        return "0/1"


def maximum_qual(records):
    return max([record.QUAL for record in records if record.QUAL is not None], default=None)


def maximum_tranche(records):
    tranches = set([get_tranche_2(record) for record in records])
    if "HIGH" in tranches:
        return "HIGH"
    elif "INTERMEDIATE" in tranches:
        return "INTERMEDIATE"
    elif "LOW" in tranches:
        return "LOW"
    else:
        return None


def get_sample_call(sample_name, records):
    """
    This function generates the Call for a single sample at at a given location, given a single record, multiple records or no record at all
    :param sample_name:
    :param records:
    :return:
    """
    call_data = vcfpy.OrderedDict.fromkeys(["GT", "TRANCHE2", "VAF"])

    if records:
        average_vaf = get_average_vaf(records)
        call_data["GT"] = get_gt(average_vaf)
        call_data["TRANCHE2"] = maximum_tranche(records)
        call_data["VAF"] = average_vaf

    return vcfpy.Call(sample=sample_name, data=call_data)


def get_average_vaf(records):
    return mean([float(record.INFO["BNDVAF"]) for record in records])


def generate_id(chrom, pos):
    return chrom + "_" + str(pos)
