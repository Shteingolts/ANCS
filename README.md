# Atomic Network Construction Script

## Overview
A collection of lammps-docused tools to help create, modify, read and write disordered elastic networks.


## Requirements
- Python 3.10 or newer


## Installation
1. Clone this repository:
```sh
git clone https://github.com/Shteingolts/ANCS
```

## Usage
Run the script using the following command:
```sh
python3 network.py <input_file_path> [output_file_name]
```


## Explanation
The script applies the following rules to construct the network:
1. Two atoms are considered bonded if the distance between them is less than or equal to the sum of their atomic radii.
2. Atoms with less than two bonds are iteratively deleted until all atoms in the network have at least three bonds.

By default, all atoms are assigned the same atom type of 1 and a mass of 1.0 units. You can modify these values in the script if needed.

## License

This project is licensed under the [MIT License](https://github.com/Shteingolts/L2N/blob/master/LICENSE.txt).