import vcfpy


def write_output(records, file_path, sample_name_to_header, chromosome_set):
    """
    Serialises the data into a VCF file.

    :param records: the list of records to serialise. No order is assumed.
    :param file_path: path to the output file
    :param sample_name_to_header: a map from the sample names to the headers
    :param chromosome_set: the set of chromosomes selected for analysis
    :return: nothing
    """
    assert len(sample_name_to_header) > 0, "At least one sample is required"

    # Sort the output by chromosome and position. We do a sort-in-place to optimise the memory
    sorting_function = lambda record: (record.CHROM, record.POS)
    records.sort(key=sorting_function)

    writer = vcfpy.Writer.from_path(file_path, header=get_header(sample_name_to_header, chromosome_set))

    for output_record in records:
        writer.write_record(output_record)

    writer.close()


def get_header(sample_name_to_header, chromosome_set):
    """
    Returns the header of the output VCF file
    :param sample_name_to_header: a dictionary from the sample names to the headers
    :param chromosome_set: the set of chromosomes selected for analysis
    :return: a vcfpy.Header
    """
    header = vcfpy.Header()

    header.add_line(vcfpy.HeaderLine(key="fileformat", value="VCFv4.3"))

    # CONTIG headers
    first_sample_header = next(iter(sample_name_to_header.values()))
    for input_header_line in first_sample_header.lines:
        if isinstance(input_header_line, vcfpy.ContigHeaderLine):
            if chromosome_set is None or input_header_line.mapping["ID"] in chromosome_set:
                header.add_line(input_header_line)

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
    sample_names = sample_name_to_header.keys()
    header.samples = vcfpy.SamplesInfos(sample_names)

    return header
