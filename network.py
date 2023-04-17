"""
A simple utility script, which takes lammps dump output and turns it into a lammps-readable file.
v. 0.1.0
"""

from __future__ import annotations

import os
import sys
from copy import deepcopy
from typing import List, TextIO


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

    def move(self, x: float = 0.0, y: float = 0.0, z: float = 0.0) -> Atom:
        self.x += x
        self.y += y
        self.z += z
        return self

    def translate(self, box: Box, direction: tuple = (0, 0, 0)) -> Atom:

        x_move = direction[0] * box.x
        y_move = direction[1] * box.y
        z_move = direction[2] * box.z
        return self.move(x_move, y_move, z_move)

    def dist(self, atom: Atom) -> float:
        return (
            (self.x - atom.x) ** 2 + (self.y - atom.y) ** 2 + (self.z - atom.z) ** 2
        ) ** 0.5

    def within_box(self, box: Box) -> bool:
        x_min, x_max = min(box.x1, box.x2), max(box.x1, box.x2)
        y_min, y_max = min(box.y1, box.y2), max(box.y1, box.y2)
        z_min, z_max = min(box.z1, box.z2), max(box.z1, box.z2)

        if (
            x_min <= self.x <= x_max
            and y_min <= self.y <= y_max
            and z_min <= self.z <= z_max
        ):
            return True
        else:
            return False

    def on_edge(self, box: Box, delta: float) -> bool:
        delta_x = delta_y = delta_z = delta
        if box.x1 == box.x2 == 0.0:
            delta_x = 0.0
        if box.y1 == box.y2 == 0.0:
            delta_y = 0.0
        if box.z1 == box.z2 == 0.0:
            delta_z = 0.0

        smaller_box = Box(
            box.x1 + delta_x,
            box.x2 - delta_x,
            box.y1 + delta_y,
            box.y2 - delta_y,
            box.z1 + delta_z,
            box.z2 - delta_z,
        )
        bigger_box = Box(
            box.x1 - delta_x,
            box.x2 + delta_x,
            box.y1 - delta_y,
            box.y2 + delta_y,
            box.z1 - delta_z,
            box.z2 + delta_z,
        )

        if self.within_box(bigger_box) and not self.within_box(smaller_box):
            return True
        else:
            return False


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
        if {self.atom1, self.atom2} == {other.atom1, other.atom2}:
            if round(self.length, 6) == round(other.length, 6):
                return True
        else:
            return False

    def __hash__(self) -> int:
        if self.atom1.atom_id > self.atom2.atom_id:
            return hash((self.atom2.atom_id, self.atom1.atom_id, round(self.length, 6)))
        else:
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
        atoms: List[Atom],
        bonds: List or None,
        box: Box,
        angles: List or None = None,
        dihedrals: List or None = None,
        impropers: List or None = None,
        atom_types: int = 1,
        bond_types: int = 0,
        angle_types: int = 0,
        dihedral_types: int = 0,
        improper_types: int = 0,
    ):

        self.box_dimensions = box.dimensions
        self.atoms = len(atoms)
        self.bonds = len(bonds) if (bonds is not None) else 0
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


class Box:
    x1: float
    x2: float
    y1: float
    y2: float
    z1: float
    z2: float

    def __init__(self, x1, x2, y1, y2, z1, z2):
        self.x1, self.x2 = min(x1, x2), max(x1, x2)
        self.y1, self.y2 = min(y1, y2), max(y1, y2)
        self.z1, self.z2 = min(z1, z2), max(z1, z2)

    def __repr__(self) -> str:
        return f"Box(x: ({self.x1}, {self.x2})\ny: ({self.y1}, {self.y2})\nz: ({self.z1}, {self.z2}))"

    def resize(self, delta: float) -> Box:
        self.x1 += delta / 2
        self.x2 -= delta / 2
        self.y1 += delta / 2
        self.y2 -= delta / 2
        self.z1 += delta / 2
        self.z2 -= delta / 2
        return self

    @property
    def x(self):
        return abs(self.x2 - self.x1)

    @property
    def y(self):
        return abs(self.y2 - self.y1)

    @property
    def z(self):
        return abs(self.z2 - self.z1)

    @property
    def dimensions(self):
        return (self.x1, self.x2, self.y1, self.y2, self.z1, self.z2)


def get_atoms(file_contents: List[str]) -> List[Atom]:
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
            atoms_start_line = i + 2
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


def delete_dangling(atoms: List[Atom]) -> tuple(list, int):
    new_atoms = [atom for atom in atoms if atom.n_bonds > 2]
    difference = len(atoms) - len(new_atoms)
    for atom in new_atoms:
        atom.n_bonds = 0
    return (new_atoms, difference)


def get_box(file_content: List[str]) -> Box:
    # get box size (x and y only for now) from lammps data file
    for i, line in enumerate(file_content):
        if "xlo" in line:
            x1, x2 = float(line.split()[0]), float(line.split()[1])
            y1, y2 = (
                float(file_content[i + 1].split()[0]),
                float(file_content[i + 1].split()[1]),
            )
            z1, z2 = (
                float(file_content[i + 2].split()[0]),
                float(file_content[i + 2].split()[1]),
            )
    try:
        return Box(x1, x2, y1, y2, z1, z2)
    except NameError:
        print("Failure reading box dimensions, probably not a valid input file.")
        sys.exit(1)


def add_spaces(string: str, width: int, indent: str = "right") -> str:
    """
    If the string is longer than provided width, returns the original string without change.
    """
    if width <= len(string):
        return string
    spaces_to_add = (width - len(string)) * " "
    if indent == "right":
        return spaces_to_add + string
    if indent == "left":
        return string + spaces_to_add


def table_row(items: List, widths: List, indent: str = "right") -> str:
    line = []
    for item, width in zip(items, widths):
        line.append(add_spaces(str(item), width, indent))

    return "".join(line) + "\n"


def make_surrounding(atoms: List[Atom], box: Box, dimensions: int = 2) -> List[Atom]:
    surrounding_atoms = set()
    # spawn neighboring atoms along the x and y axis of the bounding box, 8 in total
    if dimensions == 2:
        for atom in atoms:
            for x in (-1, 0, 1):
                for y in (-1, 0, 1):
                    if not x * box.x == y * box.y == 0:
                        surrounding_atoms.add(deepcopy(atom).translate(box, (x, y, 0)))
    elif dimensions == 3:
        for atom in atoms:
            for x in (-1, 0, 1):
                for y in (-1, 0, 1):
                    for z in (-1, 0, 1):
                        if not x * box.x == y * box.y == z * box.z == 0:
                            surrounding_atoms.add(
                                deepcopy(atom).translate(box, (x, y, z))
                            )
    else:
        raise NotImplementedError(
            "Allowed values for `dimensions` argument are either 2 or 3."
        )
    return list(surrounding_atoms)


def _make_inner_bonds(atoms: List[Atom]) -> List[Bond]:
    bonds = set()
    for atom_k in atoms:
        for atom_j in atoms:
            if atom_k != atom_j:
                if atom_k.dist(atom_j) <= (
                    (atom_k.diameter / 2) + (atom_j.diameter / 2)
                ):
                    bonds.add(Bond(atom_k, atom_j))
                    atom_k.n_bonds += 1
                    atom_j.n_bonds += 1
    return list(bonds)


def _make_outside_bonds(atoms: List[Atom], box: Box) -> List[Bond]:
    edges = [atom for atom in atoms if atom.on_edge(box, 1.0)]
    neighbors = make_surrounding(atoms, box)
    edge_neighbors = [atom for atom in neighbors if atom.on_edge(box, 1.0)]

    extra_bonds = set()
    for main_atom in edges:
        for outside_atom in edge_neighbors:
            if main_atom != outside_atom:
                # if the interatomic distance is less than the sum of two atomic radii,
                # the bond is formed
                if main_atom.dist(outside_atom) <= (
                    (main_atom.diameter / 2) + outside_atom.diameter / 2
                ):
                    extra_bonds.add(Bond(main_atom, outside_atom))
                    main_atom.n_bonds += 1
    return list(extra_bonds)


def make_bonds(atoms, box) -> list:
    bonds = set()
    for atom_k in atoms:
        for atom_j in atoms:
            if atom_k != atom_j:
                if atom_k.dist(atom_j) <= (
                    (atom_k.diameter / 2) + (atom_j.diameter / 2)
                ):
                    bonds.add(Bond(atom_k, atom_j))
                    atom_k.n_bonds += 1

    edges = [atom for atom in atoms if atom.on_edge(box, 1.0)]
    neighbors = make_surrounding(atoms, box)
    edge_neighbors = [atom for atom in neighbors if atom.on_edge(box, 1.0)]

    extra_bonds = set()
    for main_atom in edges:
        for outside_atom in edge_neighbors:
            if main_atom.dist(outside_atom) <= (
                (main_atom.diameter / 2) + outside_atom.diameter / 2
            ):
                extra_bonds.add(Bond(main_atom, outside_atom))
                main_atom.n_bonds += 1
    return list(bonds.union(extra_bonds))


def write_atoms(file: TextIO, atoms: List[Atom]):
    # 7-7-7-11-11-11-11
    legend = ["atomID", "atomType", "diameter", "density", "x", "y", "z"]
    file.write(f"\nAtoms # {legend}\n\n")
    for atom in atoms:
        properties = [
            atom.atom_id,
            "1",
            "1",
            "0.000000",
            round(atom.x, 6),
            round(atom.y, 6),
            round(atom.z, 6),
        ]
        widths = [7, 7, 7, 11, 11, 11, 11]
        line = table_row(properties, widths)
        file.write(line)


def write_bonds(file: TextIO, bonds: List[Bond]):
    # 10-10-10-10
    legend = ["ID", "type", "atom1", "atom2"]
    file.write(f"\nBonds # {legend}\n\n")
    for _id, bond in enumerate(bonds):
        properties = [_id + 1, _id + 1, bond.atom1.atom_id, bond.atom2.atom_id]
        widths = [10, 10, 10, 10]
        line = table_row(properties, widths)
        file.write(line)


def write_bond_coeffs(file: TextIO, bonds: List[Bond]):
    # 7-11-11
    legend = ["bondID", "bondCoeff", "d"]
    file.write(f"\nBond Coeffs # {legend}\n\n")
    for n, bond in enumerate(bonds):
        properties = [n + 1, round(bond.bond_coefficient, 6), round(bond.length, 6)]
        widths = [7, 11, 11]
        line = table_row(properties, widths)
        file.write(line)


def main():
    usage_info = "\n[USAGE]:\n\n    ./network.py target_file [OPTIONAL] out_file.\n"
    if len(sys.argv) < 2:
        print("\n[ERROR]: target file was not provided.")
        print(usage_info)
        sys.exit(0)
    elif sys.argv[1] == "help":
        print(usage_info)
        sys.exit(0)

    input_file = sys.argv[1]
    out_file = sys.argv[2] if len(sys.argv) > 2 else "output.lmp"

    with open(input_file, "r", encoding="utf8") as f:
        content = f.readlines()

    box = get_box(content)
    atoms = get_atoms(content)
    bonds = make_bonds(atoms, box)

    # we assume that there's at least one dangling bead
    # if not, nothing bad will happen anyway
    dangling_beads: int = 1
    steps: int = 1
    while dangling_beads > 0:
        atoms, dangling_beads = delete_dangling(atoms)
        bonds = make_bonds(atoms, box)
        steps += 1

    print(f"Atoms: {len(atoms)}.")
    print(f"Bonds: {len(bonds)}.")

    header = Header(atoms, bonds, box)

    with open(out_file, "w", encoding="utf8") as f:
        header.write_header(f)
        write_atoms(f, atoms)
        write_bonds(f, bonds)
        write_bond_coeffs(f, bonds)

    if os.path.isfile(out_file) and os.stat(out_file).st_size != 0:
        print(f"Output was written in: {out_file}")


main()
