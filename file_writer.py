import vcfpy


def write_output(output_list, output_file_path, sample_names_to_header):
    """
    Serialises the data into a VCF file.

    :param output_list:
    :param output_file_path:
    :param sample_names_to_header:
    :return:
    """
    assert len(sample_names_to_header) > 0, "At least one sample is required"

    # Sort the output by chromosome and position
    sorting_function = lambda record: (record.CHROM, record.POS)
    output_list.sort(key=sorting_function)

    header = get_header(sample_names_to_header.keys())

    writer = vcfpy.Writer.from_path(output_file_path, header)

    for output_record in output_list:
        writer.write_record(output_record)

    writer.close()


def get_header(sample_names):
    header = vcfpy.Header()

    # INFO fields
    header.add_info_line(vcfpy.OrderedDict(
        ID="END",
        Number=1,
        Type="Integer",
        Description="Stop position of the interval"))
    header.add_info_line(vcfpy.OrderedDict(
        ID="SVTYPE",
        Number=1,
        Type="String",
        Description="Type of structural variant"))

    # FORMAT fields
    header.add_format_line(vcfpy.OrderedDict(
        ID="GT",
        Number=1,
        Type="String",
        Description="Genotype"))
    header.add_format_line(vcfpy.OrderedDict(
        ID="TRANCHE2",
        Number=1,
        Type="String",
        Description="Quality category of GRIDSS structural variant calls determined using FILTER,SRQ,AS,RAS. Values are LOW INTERMEDIATE HIGH"))
    header.add_format_line(vcfpy.OrderedDict(
        ID="VAF",
        Number=1,
        Type="Float",
        Description="VAF of this SV call, derived from BNDVAF values of BND calls used to call this SV"))

    # Samples
    header.samples = vcfpy.SamplesInfos(sample_names)

    return header
