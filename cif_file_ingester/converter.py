import sys
from pypif import pif
from .parse_cif_pmg import parse_cif

def convert(files=[], **kwargs):
    """
    Convert files into a pif
    :param files: to convert
    :param kwargs: any other arguments
    :return: the pif produced by this conversion
    """
    print('Converting {} CIFs'.format(len(files)))
    systems = []
    for f in files:
        converted_pif = parse_cif(f)
        if converted_pif:
            systems.append(converted_pif)

    return systems


if __name__ == '__main__':
    with open(sys.argv[1].replace('.cif','-pif.json'), 'w') as output_file:
        pif.dump(convert(files=[sys.argv[1]]), output_file, indent=4)
