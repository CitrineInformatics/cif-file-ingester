from pypif.obj import *
import re
import ase.io
from pymatgen import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from io import open

def get_crystal_system(structure):

    '''Converts space group number to a crystal system

        Args:
            structure (structure): Pymatgen structure object

        Returns:
            sgroup (str): space group
            num (int): space group number
            bravais (str): crystal system

    '''

    sgroup, num = structure.get_space_group_info()
    if num <= 2:
        bravais = 'Triclinic'
    elif num <= 15:
        bravais = 'Monoclinic'
    elif num <= 74:
        bravais = 'Orthorhombic'
    elif num <= 142:
        bravais = 'Tetragonal'
    elif num <= 167:
        bravais = 'Trigonal'
    elif num <= 194:
        bravais = 'Hexagonal'
    else:
        bravais = 'Cubic'
    return sgroup, num, bravais

def parse_text(f):
    system = ChemicalSystem()
    system.names = []
    system.properties = []
    system.references = []
    conditions_for_some_props = []
    props = ["Lattice parameter (a)", "Lattice parameter (b)", "Lattice parameter (c)", "Cell angle ($\\alpha$)",
             "Cell angle ($\\beta$)", "Cell angle ($\\gamma$)", "Unit cell volume", "Density", "Crystal color"]
    ##### Scan through the CIF to get remaining information #####
    if f.endswith('.cif') or f.endswith('.mcif'):
        with open(f, 'r',encoding="ISO-8859-1") as cif_file:
            try:
                lines = cif_file.readlines()
            except UnicodeDecodeError:
                return system
    else:
        raise IOError('Filetype not compatible with parser. Please upload a .cif or .mcif file.\n')
    
    k_maximal = False
    for line in lines:
        if '_chemical_formula_sum' in line and '?' not in line:
            system.chemical_formula = re.sub(r'\s+','',line.split('_chemical_formula_sum')[1].strip()).strip("'")
        if "_chemical_name_systematic" in line and '?' not in line:
            system.names.append(line.split("_chemical_name_systematic")[1].strip().replace("'", ""))
        if "_journal_paper_doi" in line:
            system.references.append(Reference(doi=line.split("_journal_paper_doi")[1].strip()))
        if "_citation_DOI" in line and '?' not in line:
            system.references.append(Reference(doi=line.split("_citation_DOI")[1].strip()))
        if "_exptl_crystal_colour" in line:
            system.properties.append(Property(name="Crystal color", scalars=[Scalar(line.split("_exptl_crystal_colour")[1].strip().lower().replace('colourless', 'colorless').replace("'", ""))]))
        if '_transition_temperature' in line:
            system.properties.append(Property(name='Transition Temperature',scalars=[Scalar(value=line.split('_transition_temperature')[1].strip())],units='K'))
        if ('k-maximal' in line or 'k maximal' in line) and (f.endswith('.mcif')):
            k_maximal = True
            system.properties.append(Property(name='k-maximal subgroup - classifier',scalars=[Scalar(value=1)]))
        if '_space_group.magn_ssg_name' in line:
            system.properties.append(Property(name='Magnetic Superspace Group',scalars=[Scalar(value=line.split('_space_group.magn_ssg_name')[1].strip().strip('"'))]))
        if '_space_group.magn_point_group_name' in line or '_space_group_magn.point_group_name' in line :
            system.properties.append(Property(name='Magnetic Point Group',scalars=[Scalar(value=line.split('point_group_name')[1].strip().strip('"'))]))
        if '_space_group.magn_point_group ' in line:
            system.properties.append(Property(name='Magnetic Point Group',scalars=[Scalar(value=line.split('_space_group.magn_point_group ')[1].strip().strip('"'))]))
        if '_space_group_magn.name_BNS' in line or '_space_group.magn_name_BNS' in line:
            system.properties.append(Property(name='BNS Magnetic Point Group',scalars=[Scalar(value=re.sub(r'\s+','',line.split('name_BNS')[1]).strip('"'))]))
        if '_transition_temperature' in line or '_Neel_temperature' in line:
            system.properties.append(Property(name='Transition Temperature',scalars=[Scalar(value=line.split('_temperature')[1].strip())],units='K'))

        # Get conditions for previously specified properties
        if "_cell_measurement_temperature" in line:
            conditions_for_some_props.append(Value(name="Temperature", scalars=[Scalar(line.split("_cell_measurement_temperature")[1].strip())], units="K"))
        if '_experiment_temperature' in line:
            conditions_for_some_props.append(Value(name='Experiment Temperature',scalars=[Scalar(value=line.split('_experiment_temperature')[1].strip())],units='K'))
        if "_diffrn_radiation_type" in line:
            if "MoK\\a" in line.split("_diffrn_radiation_type")[1].strip():
                conditions_for_some_props.append(Value(name="Radiation type", scalars=[Scalar("Mo K-$\\alpha$")]))
            elif "CuK\\a" in line.split("_diffrn_radiation_type")[1].replace(" ", ""):
                conditions_for_some_props.append(Value(name="Radiation type", scalars=[Scalar("Cu K-$\\alpha$")]))
            else:
                conditions_for_some_props.append(Value(name="Radiation type", scalars=[Scalar(line.split("_diffrn_radiation_type")[1].strip())]))
        if "_diffrn_radiation_wavelength" in line:
            if line.split("_diffrn_radiation_wavelength")[1].strip().startswith("."):
                conditions_for_some_props.append(Value(name="Radiation wavelength", scalars=[Scalar("0"+line.split("_diffrn_radiation_wavelength")[1].strip())], units="$\AA$"))
            else:
                conditions_for_some_props.append(Value(name="Radiation wavelength", scalars=[Scalar(line.split("_diffrn_radiation_wavelength")[1].strip())], units="$\AA$"))
    if not k_maximal and f.split('.')[-1] == 'mcif':
        system.properties.append(Property(name='k-maximal subgroup - classifier',scalars=[Scalar(value=0)]))

    system = parse_with_pmg(f,system,lines,conditions_for_some_props)

    if len(conditions_for_some_props):
        for prop in system.properties:
            if prop.name in props:
                prop.conditions = conditions_for_some_props

    return system

def parse_with_pmg(f,system,lines,conditions_for_some_props):

    '''Loads cif file as a pymatgen structure to obtain some properties. Also reads lines in input file for properties
    that are not generated from pymatgen.

    Args:
        f (str): The input file (.cif)


    Returns:
        system: ChemicalSystem() with properties extracted from .cif file.
    '''

    tol = 3

    ##### Use pymatgen to quickly obtain some information #####
    try:
        structure = Structure.from_str(''.join(lines), fmt='cif')
    except Exception:
        # Try with ASE instead
        try:
            ase_res = ase.io.read(f,format='cif')
            if not ase_res:
                raise ValueError("ASE retrieved no data from file")
            # Convert ASE Atoms to Pymatgen Structure
            structure = AseAtomsAdaptor.get_structure(ase_res)
        except Exception:
            print(f, 'is not interpretable by pymatgen or ase')
            return system

    sg, num, bravais = get_crystal_system(structure)
    a, b, c = structure.lattice.abc
    alpha, beta, gamma = structure.lattice.angles

    system.chemical_formula = structure.composition.reduced_formula



    system.properties.append(Property(name='Space group symbol', scalars=[Scalar(value=sg)]))
    system.properties.append(Property(name='Space group number', scalars=[Scalar(value=num)]))
    system.properties.append(Property(name='Crystal system', scalars=[Scalar(value=bravais)]))
    system.properties.append(Property(name='Lattice parameter (a)', scalars=[Scalar(value=round(a, tol))], units="$\AA$"))
    system.properties.append(Property(name='Lattice parameter (b)', scalars=[Scalar(value=round(b, tol))], units="$\AA$"))
    system.properties.append(Property(name='Lattice parameter (c)', scalars=[Scalar(value=round(c, tol))], units="$\AA$"))
    system.properties.append(Property(name='Unit cell volume', scalars=[Scalar(value=round(structure.lattice.volume, tol))], units="$\AA^3$"))
    system.properties.append(Property(name='Cell angle ($\\alpha$)', scalars=[Scalar(value=round(alpha, tol))], units="$^\\circ$"))
    system.properties.append(Property(name='Cell angle ($\\beta$)', scalars=[Scalar(value=round(beta, tol))], units="$^\\circ$"))
    system.properties.append(Property(name='Cell angle ($\\gamma$)', scalars=[Scalar(value=round(gamma, tol))], units="$^\\circ$"))
    system.properties.append(Property(name='Number of atoms in unit cell', scalars=[Scalar(value=structure.composition.num_atoms)]))
    system.properties.append(Property(name='Charge of unit cell', scalars=[Scalar(value=structure.charge)]))
    system.properties.append(Property(name='Density', scalars=[Scalar(value=round(structure.density, tol))], units="g/cm$^3$"))
    system.properties.append(Property(name='Input file', files=[FileReference(relative_path=f)]))

    return system

