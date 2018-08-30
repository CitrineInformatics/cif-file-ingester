from setuptools import setup, find_packages

setup(name='cif_file_ingester',
    version='1.2.0',
    url='http://github.com/CitrineInformatics/cif_file_ingester',
    description='Converts standardized crystallographic information files (.cif) as popularized by the International Union of Crystallography (IUoC).',
    author='Chris Borg, Enze Chen, Max Hutchinson',
    author_email='cborg@citrine.io',
    packages=find_packages(),
    install_requires=[
        'pypif>=2.1.0,<3',
        'ase',
        'pymatgen'
    ],
    entry_points={
        'citrine.dice.converter': [
            'cif = cif_file_ingester.converter',
        ],
    },
)
