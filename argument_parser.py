import argparse

def parse_arguments(args):
    parser = argparse.ArgumentParser(description='Process SVs and BNDs in Gridds Format')
    parser.add_argument(dest='files', nargs='*',
                        help='Space separated files')
    parser.add_argument('-o', dest='outputFile',required=True,
                        help='Path for the output file')
    return parser.parse_args()