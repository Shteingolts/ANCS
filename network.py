"""
A simple utility script, which takes lammps dump output and turns it
into a lammps-readable file.
v. 0.1.0
"""

from __future__ import annotations

import os
import sys
from copy import deepcopy
from math import acos, degrees, sqrt
from typing import List, TextIO

from helpers import add_spaces, table_row


class Atom:
    atom_id: int
    atom_type: int
    diameter: float
    n_bonds: int
    bonded: list[int]
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
        atom_type: int = 1,
    ):
        self.atom_id = int(atom_id)
        self.atom_type = atom_type
        self.diameter = float(diameter)
        self.n_bonds = 0
        self.bonded = []
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def __repr__(self) -> str:
        return f"Atom {self.atom_id} : {self.x}, {self.y}, {self.z})."

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
        self.bond_coefficient = 1 / (self.length**2)

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


class Angle:
    """
    Value in degrees.
    """

    angle_id: int
    atom1: Atom
    atom2: Atom
    atom3: Atom
    energy: float
    value: float

    def __init__(
        self,
        angle_id: int,
        atom1: Atom,
        atom2: Atom,
        atom3: Atom,
        box: Box,
        energy: float = 0.0,
        value: float = None,
    ):
        self.angle_id = angle_id
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.energy = energy

        if value is None:
            # stupid algorithm, but works.
            # find all possible copies,
            first_atom_candidates = []
            for x in (-1, 0, 1):
                for y in (-1, 0, 1):
                    for z in (-1, 0, 1):
                        first_atom_candidates.append(
                            deepcopy(atom1).translate(box, (x, y, z))
                        )

            third_atom_candidates = []
            for x in (-1, 0, 1):
                for y in (-1, 0, 1):
                    for z in (-1, 0, 1):
                        third_atom_candidates.append(
                            deepcopy(atom3).translate(box, (x, y, z))
                        )

            def closest(origin_atom: Atom, candidates: list[Atom]):
                closest_atom = candidates[0]

                for atom in candidates[1:]:
                    if origin_atom.dist(atom) < origin_atom.dist(closest_atom):
                        closest_atom = atom
                return closest_atom

            atom1 = closest(atom2, first_atom_candidates)
            atom3 = closest(atom2, third_atom_candidates)

            v12 = [atom1.x - atom2.x, atom1.y - atom2.y, atom1.z - atom2.z]
            v23 = [atom3.x - atom2.x, atom3.y - atom2.y, atom3.z - atom2.z]

            dot_product = v12[0] * v23[0] + v12[1] * v23[1] + v12[2] * v23[2]

            mag_v12 = sqrt(v12[0] ** 2 + v12[1] ** 2 + v12[2] ** 2)
            mag_v23 = sqrt(v23[0] ** 2 + v23[1] ** 2 + v23[2] ** 2)

            cos_angle = dot_product / (mag_v12 * mag_v23)
            angle = round(degrees(acos(cos_angle)), 6)
            self.value = angle
        else:
            self.value = value

    def __eq__(self, other: Angle) -> bool:
        if self.value == other.value:
            return True
        return False

    def __hash__(self) -> int:
        if self.atom1.atom_id > self.atom3.atom_id:
            return hash((self.atom3, self.atom2, self.atom1))
        else:
            return hash((self.atom1, self.atom2, self.atom3))

    def __repr__(self) -> str:
        return (
            f"Angle {self.angle_id} : "
            f"{self.atom1.atom_id}-{self.atom2.atom_id}"
            f"-{self.atom3.atom_id} | {round(self.value, 2)} "
            f"({round(180 - self.value, 2)}) deg."
        )


class Dihedral:
    atom1: Atom
    atom2: Atom
    atom3: Atom
    atom4: Atom

    def __init__(self) -> None:
        raise NotImplementedError("not yet...")


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
        atom_types: int = 1,  # defaults to one
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
        return (
            f"Box ({round(self.x1, 3)} : {round(self.x2, 3)}) "
            f"({round(self.y1, 3)} : {round(self.y2, 3)}) "
            f"({round(self.z1, 3)} : {round(self.z2, 3)})"
        )

    @classmethod
    def from_atoms(cls, file: str) -> Box:
        with open(file, encoding="utf8") as atoms_file:
            content = atoms_file.readlines()
            x1, x2 = (float(content[5].split()[0]), float(content[5].split()[1]))
            y1, y2 = (
                float(content[6].split()[0]),
                float(content[6].split()[1]),
            )
            z1, z2 = (
                float(content[7].split()[0]),
                float(content[7].split()[1]),
            )
        try:
            return Box(x1, x2, y1, y2, z1, z2)
        except NameError:
            print("Failure reading box dimensions, probably not a valid input file.")
            sys.exit(1)

    @classmethod
    def from_data_file(cls, file: str) -> Box:
        with open(file, encoding="utf8") as data_file:
            content = data_file.readlines()

        x1, x2 = (
            float(content[11].split()[0]),
            float(content[11].split()[1]),
        )
        y1, y2 = (
            float(content[12].split()[0]),
            float(content[12].split()[1]),
        )
        z1, z2 = (
            float(content[13].split()[0]),
            float(content[13].split()[1]),
        )
        return Box(x1, x2, y1, y2, z1, z2)

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


class Network:
    atoms: list[Atom]
    bonds: list[Bond]
    angles: list[Angle] or None
    dihedrals: list[Dihedral] or None
    masses: dict[int:float] or None
    box: Box
    header: Header

    def __init__(
        self,
        atoms: list[Atom],
        bonds: list[Bond],
        box: Box,
        header: Header,
        angles: list[Angle] = None,
        dihedrals: list[Dihedral] = None,
        masses: dict or None = None,
    ):
        self.atoms = atoms
        self.bonds = bonds
        self.box = box
        self.header = header
        self.angles = angles
        self.dihedrals = dihedrals
        self.masses = masses

    @staticmethod
    def _compute_angles(atoms: list[Atom], box: Box) -> list[Angle]:
        atoms_map = {atom.atom_id: atom for atom in atoms}
        angles = set()
        for atom in atoms:
            for neighbour_k in atom.bonded:
                for neighbour_j in atom.bonded:
                    if atoms_map[neighbour_k] != atoms_map[neighbour_j]:
                        angles.add(
                            Angle(
                                len(angles) + 1,
                                atoms_map[neighbour_k],
                                atom,
                                atoms_map[neighbour_j],
                                box,
                            )
                        )
        return list(angles)

    def _compute_dihedrals(self):
        raise NotImplementedError("not yet...")

    def remove_bond(self, bond: Bond):
        try:
            self.bonds.remove(bond)
        except ValueError:
            print(f"[ERROR] Bond {bond} is not in the network!")
            return

    @classmethod
    def from_data_file(cls, file: str) -> Network:
        """
        LAMMPS data file parser, returns a network.
        At least `Atoms` and `Bonds` sections should be present in the file.
        """
        with open(file, encoding="utf8") as data_file:
            content = data_file.readlines()

        box = Box(
            float(content[11].split()[0]),
            float(content[11].split()[1]),
            float(content[12].split()[0]),
            float(content[12].split()[1]),
            float(content[13].split()[0]),
            float(content[13].split()[1]),
        )

        n_atoms = int(content[1].split()[0])
        n_bonds = int(content[2].split()[0])
        n_angles = int(content[3].split()[0])
        n_dihedrals = int(content[4].split()[0])

        atoms = []
        bonds = []
        angles = []
        # dihedrals = []

        location = {
            "atoms": (),
            "bonds": (),
            "angles": (),
            "dihedrals": (),
        }

        for index, line in enumerate(content):
            if "Atoms" in line.strip():
                atoms_start = index + 2
                atoms_end = atoms_start + n_atoms
                location["atoms"] = (atoms_start, atoms_end)
            if "Bonds" in line.strip():
                bonds_start = index + 2
                bonds_end = bonds_start + n_bonds
                location["bonds"] = (bonds_start, bonds_end)
            if "Angles" in line.strip():
                angles_start = index + 2
                angles_end = angles_start + n_angles
                location["angles"] = (angles_start, angles_end)
            if "Dihedrals" in line.strip():
                dihedrals_start = index + 2
                dihedrals_end = dihedrals_start + n_dihedrals
                location["dihedrals"] = (dihedrals_start, dihedrals_end)

        if location["atoms"]:
            print(f"Atoms expected: {n_atoms}")
            for line in content[atoms_start:atoms_end]:
                data = line.split()
                atoms.append(
                    Atom(
                        int(data[0]),
                        0.0,
                        float(data[4]),
                        float(data[5]),
                        float(data[6]),
                    )
                )
            print(f"Atoms read: {len(atoms)}")
        else:
            print("[ERROR]: Something went wrong when reading atoms from the file.")

        if location["bonds"]:
            print(f"Bonds expected: {n_bonds}")
            atoms_map = {atom.atom_id: atom for atom in atoms}
            for line in content[bonds_start:bonds_end]:
                data = line.split()
                atom1_id = int(data[2])
                atom2_id = int(data[3])
                atom1 = atoms_map[atom1_id]
                atom2 = atoms_map[atom2_id]
                atom1.bonded.append(atom2_id)
                atom2.bonded.append(atom1_id)
                bonds.append(Bond(atom1, atom2))

            print(f"Bonds read: {len(bonds)}")
        else:
            print("[ERROR]: Something went wrong when reading bonds from the file.")

        # at this point, the bare minumum for the network sould be present
        header = Header(atoms, bonds, box)
        local_network = Network(atoms, bonds, box, header)

        if location["angles"]:
            print(f"Angles expected: {n_angles}")
            atoms_map = {atom.atom_id: atom for atom in atoms}
            for line in content[angles_start:angles_end]:
                data = line.split()
                angle_id = int(data[0])
                atom1 = atoms_map[int(data[2])]
                atom2 = atoms_map[int(data[3])]
                atom3 = atoms_map[int(data[4])]
                angles.append(Angle(angle_id, atom1, atom2, atom3, box))
            print(f"Angles read: {len(angles)}")
            local_network.angles = angles
            header.angles = len(angles)
            header.angle_types = len(angles)
        else:
            print("No angle data have been found")
            print("Calculating angles..")
            angles = Network._compute_angles(atoms, box)
            local_network.angles = angles
            header.angles = len(angles)
            header.angle_types = len(angles)
            print(f"Angles calculated: {len(angles)}")

        if location["dihedrals"]:
            print(f"Dihedrals expected: {n_dihedrals}")
            print("Dihedrals are not yet emplemented.")
            # TODO Implement reading and writing dihedrals
        else:
            print("No dihedrals data have been found")

        return local_network

    @classmethod
    def from_atoms(cls, input_file: str) -> Network:
        """
        Computes bonds from the atomic coordinates.
        """
        with open(input_file, "r", encoding="utf8") as f:
            content = f.readlines()

        box = Box.from_atoms(input_file)
        print(f"{box}")
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

        print(f"Atoms: {len(atoms)}")
        print(f"Bonds: {len(bonds)}")

        angles = Network._compute_angles(atoms, box)
        header = Header(atoms, bonds, box, angles=angles)
        return Network(atoms, bonds, box, header, angles=angles)

    def write_to_file(self, target_file: str) -> str:
        """
        Writes network to a file a returns the path
        """
        path = os.path.abspath(os.path.join(os.getcwd(), target_file))
        with open(path, "w", encoding="utf8") as file:
            self.header.write_header(file)

            if self.atoms:
                # write `Atoms` section
                # 7-7-7-11-11-11-11
                legend = ["atomID", "moleculeID", "atomType", "charge", "x", "y", "z"]
                file.write(f"\nAtoms # {legend}\n\n")
                for atom in self.atoms:
                    properties = [
                        atom.atom_id,
                        "1",  # always 1 for now
                        atom.atom_type,  # defaults to 1 when construsted
                        "0.000000",  # always neutral for now
                        round(atom.x, 6),
                        round(atom.y, 6),
                        round(atom.z, 6),
                    ]
                    widths = [7, 7, 7, 11, 11, 11, 11]
                    line = table_row(properties, widths)
                    file.write(line)

            file.write("\n# Masses\n\n")
            file.write("# Uncomment and replace this with proper values by hand\n")

            if self.bonds:
                # write `Bonds` section
                # 10-10-10-10
                legend = ["ID", "type", "atom1", "atom2"]
                file.write(f"\nBonds # {legend}\n\n")
                for _id, bond in enumerate(self.bonds):
                    properties = [
                        _id + 1,
                        _id + 1,
                        bond.atom1.atom_id,
                        bond.atom2.atom_id,
                    ]
                    widths = [10, 10, 10, 10]
                    line = table_row(properties, widths)
                    file.write(line)

                # write `Bond Coeffs` section
                # 7-11-11
                legend = ["bondID", "bondCoeff", "d"]
                file.write(f"\nBond Coeffs # {legend}\n\n")
                for n, bond in enumerate(self.bonds):
                    properties = [
                        n + 1,
                        round(bond.bond_coefficient, 6),
                        round(bond.length, 6),
                    ]
                    widths = [7, 11, 11]
                    line = table_row(properties, widths)
                    file.write(line)

            if self.angles:
                # write `Angles` section
                # 10-10-10-10-10
                legend = ["angleID", "angleType", "atom1", "atom2", "atom3"]
                file.write(f"\nAngles # {legend}\n\n")
                for angle in self.angles:
                    properties = [
                        angle.angle_id,
                        angle.angle_id,  # type the same as id
                        angle.atom1.atom_id,
                        angle.atom2.atom_id,
                        angle.atom3.atom_id,
                    ]
                    widths = [10, 10, 10, 10, 10]
                    line = table_row(properties, widths)
                    file.write(line)

                # write `Angle Coeffs` section
                # 7-11-11
                legend = ["angleID", "energy", "value (deg)"]
                file.write(f"Angle Coeffs # {legend}\n\n")
                for angle in self.angles:
                    properties = [angle.angle_id, angle.energy, angle.value]
                    widths = [7, 11, 11]
                    line = table_row(properties, widths)
                    file.write(line)

            if self.dihedrals:
                print("Writing dihedrals to a file is not implemented.")
                # TODO Implement reading and writing dihedrals

        if os.path.exists(path) and os.path.getsize(path) > 0:
            print(f"Output was written in: {os.path.abspath(path)}")
            return path
        print(f"Problem saving network to {path}")
        sys.exit(1)


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


def make_bonds(atoms: list[Atom], box: Box) -> list:
    bonds = set()
    for atom_k in atoms:
        for atom_j in atoms:
            if atom_k != atom_j:
                if atom_k.dist(atom_j) <= (
                    (atom_k.diameter / 2) + (atom_j.diameter / 2)
                ):
                    bonds.add(Bond(atom_k, atom_j))
                    atom_k.bonded.append(atom_j.atom_id)
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
                main_atom.bonded.append(outside_atom.atom_id)
                main_atom.n_bonds += 1
    return list(bonds.union(extra_bonds))


def main():
    usage_info = "\n[USAGE]:\n\n    ./network.py target_file [OPTIONAL] out_file.\n"
    if len(sys.argv) < 2:
        print("\n[ERROR]: target file was not provided.")
        print(usage_info)
        sys.exit(0)
    elif sys.argv[1] == "help":
        print(usage_info)
        sys.exit(0)

    input_file_path = sys.argv[1]
    print(f"Input file: {os.path.abspath(input_file_path)}")

    if len(sys.argv) > 2:
        out_file_path = sys.argv[2]
    else:
        input_file_path = os.path.abspath(input_file_path)
        input_dir = os.path.dirname(input_file_path)
        input_file_name = os.path.basename(input_file_path).split(".")[0]
        out_file_name = "".join((input_file_name, "_out.lmp"))
        out_file_path = os.path.join(input_dir, out_file_name)

    # constructing the bare minimum network from atomic coordinates
    network = Network.from_atoms(input_file_path)

    network.write_to_file(out_file_path)


main()
