from __future__ import annotations

import os
import sys
from typing import TextIO


class Atom:
    atom_id: int
    diameter: float
    n_bonds: int
    x: float
    y: float
    z: float

    def __init__(
        self,
        atom_id: int,
        diameter: float,
        x: float or int,
        y: float or int,
        z: float or int,
    ):
        self.atom_id = int(atom_id)
        self.diameter = float(diameter)
        self.n_bonds = 0
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def __repr__(self) -> str:
        return f"Atom(id: {self.atom_id}, x: {self.x}, y: {self.y}, z: {self.z})"

    def __eq__(self, other: Atom) -> bool:
        if (
            self.atom_id == other.atom_id
            and self.diameter == other.diameter
            and self.x == other.x
            and self.y == other.y
            and self.z == other.z
        ):
            return True
        return False

    def __hash__(self) -> int:
        return hash((self.atom_id, self.x, self.y, self.z))

    def dist(self, atom: Atom) -> float:
        return (
            (self.x - atom.x) ** 2 + (self.y -
                                      atom.y) ** 2 + (self.z - atom.z) ** 2
        ) ** 0.5


class Bond:
    """
    All bonds haver the same elasticity. The bond coefficient is 1 / d^2,
    where d is the bond length.
    """

    atom1: Atom
    atom2: Atom
    length: float
    bond_coefficient: float

    def __init__(self, atom1: Atom, atom2: Atom):
        self.atom1 = atom1
        self.atom2 = atom2
        self.length = atom1.dist(atom2)
        self.bond_coefficient = 1 / (self.length ** 2)

    def __repr__(self) -> str:
        return f"""Bond(atom1: {self.atom1}, 
         atom2: {self.atom2}, 
         d: {self.length}, 
         coeff: {self.bond_coefficient})"""

    def __eq__(self, other: Bond) -> bool:
        if ({self.atom1, self.atom2} == {other.atom1, other.atom2}
                and round(self.length, 6) == round(other.length, 6)):
            return True
        else:
            return False

    def __hash__(self) -> int:
        return hash((self.atom1.atom_id, self.atom2.atom_id, round(self.length, 6)))


class Header:
    atoms: int
    bonds: int
    angles: int
    dihedrals: int
    impropers: int
    atom_types: int
    bond_types: int
    angle_types: int
    dihedral_types: int
    improper_types: int
    box_dimensions: tuple

    def __init__(
        self,
        atoms: list[Atom],
        bonds: list[Bond],
        box_dimensions: tuple,
        angles: list or None = None,
        dihedrals: list or None = None,
        impropers: list or None = None,
        atom_types: int = 1,
        bond_types: int = 0,
        angle_types: int = 0,
        dihedral_types: int = 0,
        improper_types: int = 0,
    ):

        self.atoms = len(atoms)
        self.bonds = len(bonds)
        self.angles = len(angles) if (angles is not None) else 0
        self.dihedrals = len(dihedrals) if (dihedrals is not None) else 0
        self.impropers = len(impropers) if (impropers is not None) else 0
        self.atom_types = atom_types

        if not self.bonds:
            self.bond_types = 0
        else:
            if bond_types == 0:
                self.bond_types = len(bonds)
            else:
                self.bond_types = bond_types

        if not self.angles:
            self.angle_types = 0
        else:
            if angle_types == 0:
                self.angle_types = len(angles)
            else:
                self.angle_types = angle_types

        if not self.dihedrals:
            self.dihedral_types = 0
        else:
            if dihedral_types == 0:
                self.dihedral_types = len(dihedrals)
            else:
                self.dihedral_types = dihedral_types

        if not self.impropers:
            self.improper_types = 0
        else:
            if improper_types == 0:
                self.improper_types = len(improper_types)
            else:
                self.improper_types = improper_types

        self.box_dimensions = box_dimensions

    def write_header(self, file: TextIO) -> None:
        file.write("LAMMPS data file.\n")
        file.write(add_spaces(f"{str(self.atoms)}", 7))
        file.write(" atoms\n")
        file.write(add_spaces(f"{str(self.bonds)}", 7))
        file.write(" bonds\n")
        file.write(add_spaces(f"{str(self.angles)}", 7))
        file.write(" angles\n")
        file.write(add_spaces(f"{str(self.dihedrals)}", 7))
        file.write(" dihedrals\n")
        file.write(add_spaces(f"{str(self.impropers)}", 7))
        file.write(" impropers\n")
        file.write(add_spaces(f"{str(self.atom_types)}", 7))
        file.write(" atom types\n")
        file.write(add_spaces(f"{str(self.bond_types)}", 7))
        file.write(" bond types\n")
        file.write(add_spaces(f"{str(self.angle_types)}", 7))
        file.write(" angle types\n")
        file.write(add_spaces(f"{str(self.dihedral_types)}", 7))
        file.write(" dihedral types\n")
        file.write(add_spaces(f"{str(self.improper_types)}", 7))
        file.write(" improper types\n")
        file.write(add_spaces(f"{str(round(self.box_dimensions[0], 6))}", 11))
        file.write(add_spaces(f"{str(round(self.box_dimensions[1], 6))}", 11))
        file.write(" xlo xhi\n")
        file.write(add_spaces(f"{str(round(self.box_dimensions[2], 6))}", 11))
        file.write(add_spaces(f"{str(round(self.box_dimensions[3], 6))}", 11))
        file.write(" ylo yhi\n")
        file.write(add_spaces(f"{str(round(self.box_dimensions[4], 6))}", 11))
        file.write(add_spaces(f"{str(round(self.box_dimensions[5], 6))}", 11))
        file.write(" zlo zhi\n")


def get_atoms(file_contents: list[str]) -> list[Atom]:
    n_atoms = 0
    atoms_start_line = 0
    atoms_end_line = 0
    atoms = []

    for i, line in enumerate(file_contents[:20]):
        # skip the comment lines
        if line[0] == "#":
            continue
        # get rid the inline comments
        if "#" in line:
            line = " ".join(line.split()[: line.find("#")])
        # read number of atoms at the top of data file
        if "atoms" in line.split():
            n_atoms = int(line.split()[0])
        # find the Atoms part
        if "Atoms" in line.split():
            atoms_start_line = i + 3
            atoms_end_line = atoms_start_line + n_atoms
            break
    # Go line-by-line extracting useful info
    for atom_line in file_contents[atoms_start_line:atoms_end_line]:
        atom_id = atom_line.strip().split()[0]
        atom_diameter = atom_line.strip().split()[2]
        x = atom_line.split()[4]
        y = atom_line.split()[5]
        z = atom_line.split()[6]
        atoms.append(Atom(atom_id, atom_diameter, x, y, z))
    return atoms


def delete_dangling(atoms: list[Atom]) -> list[Atom]:
    atoms = [atom for atom in atoms if atom.n_bonds > 2]
    return atoms


def get_box(file_content: list[str]) -> tuple(float):
    # get box size (x and y only for now) from lammps data file
    for i, line in enumerate(file_content):
        if "xlo" in line:
            x1, x2 = float(line.split()[0]), float(line.split()[1])
            y1, y2 = float(file_content[i+1].split()
                           [0]), float(file_content[i+1].split()[1])
            z1, z2 = float(file_content[i+2].split()
                           [0]), float(file_content[i+2].split()[1])
            # print(f"Box: {x1}, {x2}, {y1}, {y2}, {z1}, {z2}.")
            return (x1, x2, y1, y2, z1, z2)


def add_spaces(string: str, width: int, indent: str = "right") -> str:
    if width <= len(string):
        print(f"String {string} is longer than provided width {width}!")
        return "!"
    spaces_to_add = (width - len(string)) * " "
    if indent == 'right':
        return spaces_to_add + string
    elif indent == 'left':
        return string + spaces_to_add


def table_row(items: list, widths: list, indent: str = 'right') -> str:
    line = []
    for item, width in zip(items, widths):
        line.append(add_spaces(str(item), width, indent))
    
    return ''.join(line) + '\n'


def make_bonds(atoms: list[Atom]) -> list[Bond]:
    bonds = set()
    for atom_k in atoms:
        for atom_j in atoms:
            if atom_k != atom_j:
                if atom_k.dist(atom_j) <= ((atom_k.diameter / 2) + (atom_j.diameter / 2)):  # noqa: E501
                    bonds.add(Bond(atom_k, atom_j))
                    atom_k.n_bonds += 1
                    atom_j.n_bonds += 1
    return list(bonds)


def write_atoms(file: TextIO, atoms: list[Atom]):
    # 7-7-7-11-11-11-11
    file.write("\nAtoms\n\n")
    for atom in atoms:
        properties = [
            atom.atom_id,
            '1',
            '1',
            '0.000000',
            round(atom.x, 6),
            round(atom.y, 6),
            round(atom.z, 6),
        ]
        widths = [7, 7, 7, 11, 11, 11, 11]
        line = table_row(properties, widths)
        file.write(line)


def write_bonds(file: TextIO, bonds: list[Bond]):
    # 10-10-10-10
    file.write("\nBonds\n\n")
    for _id, bond in enumerate(bonds):
        properties = [
            _id + 1,
            _id + 1,
            bond.atom1.atom_id,
            bond.atom2.atom_id
        ]
        widths = [10, 10, 10, 10]
        line = table_row(properties, widths)
        file.write(line)


def write_bond_coeffs(file: TextIO, bonds: list[Bond]):
    # 7-11-11
    file.write("\nBond Coeffs\n\n")
    for n, bond in enumerate(bonds):
        properties = [
            n+1,
            round(bond.bond_coefficient, 6),
            round(bond.length, 6)
        ]
        widths = [7, 11, 11]
        line = table_row(properties, widths)
        file.write(line)


def main():
    usage_info = """USAGE:
        ./network.py target_file [OPTIONAL] out_file."""
    if len(sys.argv) < 2:
        print("[ERROR]: target file was not provided.")
        print(usage_info)
        exit(0)
    elif sys.argv[1] == 'help':
        print(usage_info)
        exit(0)

    input_file = sys.argv[1]
    out_file = sys.argv[2] if len(sys.argv) > 2 else "output.lmp"

    with open(input_file, "r", encoding="utf8") as f:
        content = f.readlines()

    atoms = get_atoms(content)
    bonds = make_bonds(atoms)
    box = get_box(content)

    atoms = delete_dangling(atoms)
    bonds = make_bonds(atoms)
    print(f"Atoms after deletion: {len(atoms)}")
    print(f"Bonds after deletion: {len(bonds)}")

    header = Header(atoms, bonds, box)
    
    with open(out_file, "w", encoding="utf8") as f:
        header.write_header(f)
        write_atoms(f, atoms)
        write_bonds(f, bonds)
        write_bond_coeffs(f, bonds)

    if os.path.isfile(out_file) and os.stat(out_file).st_size != 0:
        print(f"Output was written in: {out_file}")


main()
