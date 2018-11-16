from cif_file_ingester.converter import convert
from cif_file_ingester.parse_cif import parse_with_pmg, get_crystal_system, parse_text
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
    system = parse_with_pmg('test_files/Al2O3.cif')
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


    try:
        parse_with_pmg('foo')
    except IOError:
        pass
    else:
        assert False, 'Bad files should return None!'

    system = parse_with_pmg('test_files/C13H22O3.cif')
    for prop in system.properties:
        if prop.conditions:
            for cond in prop.conditions:
                if cond.name == "Radiation type":
                    assert cond.scalars[0].value == 'Mo K-$\\alpha$', 'Incorrect radiation type'

    # Test mcif fields
    system = parse_with_pmg('test_files/magndata_0_53.mcif')
    for prop in system.properties:
        if prop.name == 'Magnetic Point Group':
            assert prop.scalars[0].value == "2'/m'", 'Incorrect magnetic point group'
        if prop.name == 'BNS Magnetic Point Group':
            assert prop.scalars[0].value == "C2'/m'", 'Incorrect BNS magnetic point group'
        if prop.name == 'Transition Temperature':
            assert prop.scalars[0].value == '560', 'Incorrect transition temperature'
        if prop.name == 'k-maximal subgroup - classifier':
            assert prop.scalars[0].value == 1, 'Incorrect k-maximal subgroup classifier'
    
    system = parse_with_pmg('test_files/magndata_1_1_35.mcif')
    try:
        system.properties
    except AttributeError:
        pass
    else:
        assert system is None, 'Pymatgen does not support incommensurate structures; system should be None type'

    # Test ASE fallback
    system = parse_with_pmg('test_files/Ce3VO16.cif')
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
    test_files = ['test_files/Al2O3.cif', 'test_files/C13H22O3.cif', 'test_files/C_mp-66_symmetrized.cif']
    systems = convert(test_files)
    assert len(systems) == 3, 'Incorrect number of PIFs'


def test_file_ref():
    '''
    Test that it parses CIFs correctly
    '''
    system = parse_with_pmg('test_files/C_mp-66_symmetrized.cif')
    for prop in system.properties:
        if prop.files:
            for ref in prop.files:
                assert hasattr(ref, 'relative_path') and ref.relative_path != None, "file reference not generated correctly"
        if prop.name == "Space group symbol":
            assert prop.scalars[0].value == 'Fd-3m', 'Incorrect property value'



if __name__ == '__main__':
    test_crystal_system()
    test_parse_cif()
    test_cif_converter()
