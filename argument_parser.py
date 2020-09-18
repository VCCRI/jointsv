import argparse


def parse_arguments(args):
    parser = argparse.ArgumentParser(description='Process SVs and BNDs in Gridds Format')
    parser.add_argument(dest='files', nargs='*',
                        help='Space separated files')
    parser.add_argument('-o', '--output', dest='output_file', required=False, default='output.vcf',
                        help='Path for the output file. If omitted, output.vcf will be used. If set to - then the output will be sent to stdout')
    parser.add_argument('--chromosome', dest='chromosome_list',
                        required=False, nargs='+',
                        help='Space separated chromosomes we want to analyze')
    return parser.parse_args()
