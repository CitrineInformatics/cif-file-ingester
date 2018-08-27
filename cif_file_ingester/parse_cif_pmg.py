from pypif.obj import *
from pymatgen import Structure

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


def parse_cif(f):

    '''Loads cif file as a pymatgen structure to obtain some properties. Also reads lines in input file for properties
    that are not generated from pymatgen.

    Args:
        f (str): The input file (.cif)


    Returns:
        system: ChemicalSystem() with properties extracted from .cif file.
    '''

    system = ChemicalSystem()
    system.names = []
    system.properties = []
    system.references = []

    conditions_for_some_props = []
    props = ["Lattice parameter (a)", "Lattice parameter (b)", "Lattice parameter (c)", "Cell angle ($\\alpha$)",
             "Cell angle ($\\beta$)", "Cell angle ($\\gamma$)", "Unit cell volume", "Density", "Crystal color"]
    tol = 3

    ##### Use pymatgen to quickly obtain some information #####
    try:
        structure = Structure.from_file(f)
    except:
        print(f, 'is not interpretable by pymatgen')
        return None

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


    ##### Scan through the CIF to get remaining information #####
    with open(f, 'r') as cif_file:
        lines = cif_file.readlines()
        for line in lines:
            if "_chemical_name_systematic" in line and '?' not in line:
                system.names.append(line.split("_chemical_name_systematic")[1].strip().replace("'", ""))
            if "_journal_paper_doi" in line:
                system.references.append(Reference(doi=line.split("_journal_paper_doi")[1].strip()))
            if "_exptl_crystal_colour" in line:
                system.properties.append(Property(name="Crystal color", scalars=[Scalar(line.split("_exptl_crystal_colour")[1].strip().lower().replace('colourless', 'colorless').replace("'", ""))]))

            # Get conditions for previously specified properties
            if "_cell_measurement_temperature" in line:
                conditions_for_some_props.append(Value(name="Temperature", scalars=[Scalar(line.split("_cell_measurement_temperature")[1].strip())], units="K"))
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

    if len(conditions_for_some_props):
        for prop in system.properties:
            if prop.name in props:
                prop.conditions = conditions_for_some_props

    return system

