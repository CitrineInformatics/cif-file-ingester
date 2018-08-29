from cif_file_ingester.converter import convert
from cif_file_ingester.parse_cif_pmg import parse_cif, get_crystal_system
from pypif import pif
from pypif.obj import *
from pymatgen import *
import numpy as np
import os, sys

def test_crystal_system():
    '''
    Test that the correct crystal system information is returned
    '''
    structure = Structure.from_file('test_files/Al2O3.cif')
    sg, num, bravais = get_crystal_system(structure)
    assert sg == 'R-3c', 'Incorrect space group symbol'
    assert num == 167, 'Incorrect space group number'
    assert bravais == 'Trigonal', 'Incorrect crystal system'


def test_parse_cif():
    '''
    Test that it parses CIFs correctly
    '''
    system = parse_cif('test_files/Al2O3.cif')
    assert system.chemical_formula == 'Al2O3', 'Incorrect chemical formula'
    assert system.names == ['Aluminium oxide'], 'Incorrect name'
    assert system.references[0].doi == '10.1002/pssa.2210870204'

    assert len(system.properties) == 14, 'Incorrect number of properties'
    for prop in system.properties:
        if prop.name == "Space group symbol":
            assert len(prop.scalars) == 1, 'Incorrect length'
            assert prop.scalars[0].value == 'R-3c', 'Incorrect property value'
        if prop.name == "Unit cell volume":
            assert str(prop.scalars[0].value) == '255.033', 'Incorrect unit cell volume'


    system = parse_cif('foo')
    assert system == None, 'Bad files should return None!'

    system = parse_cif('test_files/C13H22O3.cif')
    for prop in system.properties:
        if prop.conditions:
            for cond in prop.conditions:
                if cond.name == "Radiation type":
                    assert cond.scalars[0].value == 'Mo K-$\\alpha$', 'Incorrect radiation type'

    # Test ASE fallback
    system = parse_cif('test_files/Ce3VO16.cif')
    assert system.chemical_formula == 'Ce3VO16', 'Incorrect chemical formula'

    assert len(system.properties) == 14, 'Incorrect number of properties'
    for prop in system.properties:
        if prop.name == "Space group symbol":
            assert len(prop.scalars) == 1, 'Incorrect length'
            assert prop.scalars[0].value == 'I4_1/amd', 'Incorrect property value'
        elif prop.name == "Unit cell volume":
            assert str(prop.scalars[0].value) == '355.625', 'Incorrect unit cell volume'


def test_cif_converter():
    '''
    Test that correct number of PIFs are created
    '''
    test_files = ['test_files/Al2O3.cif', 'test_files/C13H22O3.cif', 'test_files/diamond.cif']
    systems = convert(test_files)
    assert len(systems) == 3, 'Incorrect number of PIFs'


def test_file_ref():
    '''
    Test that it parses CIFs correctly
    '''
    system = parse_cif('test_files/diamond.cif')
    for prop in system.properties:
        if prop.files:
            for ref in prop.files:
                assert hasattr(ref, 'relative_path') and ref.relative_path != None, "file reference not generated correctly"


if __name__ == '__main__':
    test_crystal_system()
    test_parse_cif()
    test_cif_converter()
