import vcfpy


def write_output(records, file_path, sample_names):
    """
    Serialises the data into a VCF file.

    :param records: the list of records to serialise. No order is assumed.
    :param file_path: path to the output file
    :param sample_names: list of the sample names
    :return: nothing
    """
    assert len(sample_names) > 0, "At least one sample is required"

    # Sort the output by chromosome and position. We do a sort-in-place to optimise the memory
    sorting_function = lambda record: (record.CHROM, record.POS)
    records.sort(key=sorting_function)

    writer = vcfpy.Writer.from_path(file_path, header=get_header(sample_names))

    for output_record in records:
        writer.write_record(output_record)

    writer.close()


def get_header(sample_names):
    """
    Returns the header of the output VCF file
    :param sample_names: iterable with the list of the sample names
    :return: a vcfpy.Header
    """
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
    header.add_info_line(vcfpy.OrderedDict(
        ID="INSSEQ",
        Number=1,
        Type="String",
        Description="Insertion sequence of structural variant, not including sequence marked as duplication"))

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
    header.add_format_line(vcfpy.OrderedDict(
        ID="INSSEQ",
        Number=1,
        Type="String",
        Description="Insertion sequence of structural variant, not including sequence marked as duplication"))

    # Samples
    header.samples = vcfpy.SamplesInfos(sample_names)

    return header
