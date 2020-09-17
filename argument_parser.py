import argparse


def parse_arguments(args):
    parser = argparse.ArgumentParser(description='Process SVs and BNDs in Gridds Format')
    parser.add_argument(dest='files', nargs='*',
                        help='Space separated files')
    parser.add_argument('-o', dest='output_file', required=False, default='output.vcf',
                        help='Path for the output file')
    parser.add_argument('--chromosome', dest='chromosome_list',
                        required=False, nargs='+',
                        help='Space separated chromosomes we want to analyze')
    return parser.parse_args()
