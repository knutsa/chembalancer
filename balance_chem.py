from sympy import Matrix
from diophantine import solve
import re
from collections import defaultdict

elems = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
         'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
         'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
         'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
         'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
         'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

def extract_molecules(molecule_pattern, side):
    for num_molecules in range(1,1000):
        side_pattern = "^\s*" + "\s*\+\s*".join([f"(?P<M{i}>{molecule_pattern})" for i in range(num_molecules)]) + "\s*$"
        found = re.search(side_pattern, side)
        if found:
            return [found.group(f"M{i}") for i in range(num_molecules)]

    raise Exception(f"No molecules found when trying up to 1000 i  {side}")

def extract_atom_data(molecule):
    atom_pattern = "(" + '|'.join([elem for elem in elems]) + ")(\d*)"
    for num_atoms in range(1, 10000):
        molecule_pattern = "^" + "".join([ f"(?P<A{i}>{atom_pattern})" for i in range(num_atoms)]) + "(?P<CHARGE>(\s\d+)?(\+|-))?$"
        found = re.search(molecule_pattern, molecule)
        if found:
            atom_data = defaultdict(int)
            atoms = [found.group(f"A{i}") for i in range(num_atoms)]
            for atom in atoms:
                element, count = re.search("^" + atom_pattern + "$", atom).groups()
                if count:
                    atom_data[element] += int(count)
                else:
                    atom_data[element] += 1
            if found.group("CHARGE"):
                charge_str = found.group("CHARGE")
                charge_str = charge_str.strip()
                abs_size = 1 if len(charge_str) == 1 else int(charge_str[:-1])
                if charge_str[-1] == '+':
                    atom_data["CHARGE"] = abs_size
                else:
                    atom_data["CHARGE"] = -abs_size
            
            return atom_data
        

def take_input()-> 'tuple[list, list]':
    """
        Asks user to input chemical reaction
        returns lists of the molecules at each side of the reaction
        A molecule is represented as a dictionary withe keys being the atoms in it, and values being the number it occurs
    """
    atom_pattern = "(" + '|'.join([elem for elem in elems]) + ")(\d*)"

    molecule_pattern = f"({atom_pattern}(\d+(\+|-))?)+((\s\d+)?(\+|-))?"
    reaction_side_pattern  = f"{molecule_pattern}(\s*\+\s*{molecule_pattern})*"
    reaction_pattern = f"^(?P<LH>{reaction_side_pattern})\s*=\s*(?P<RH>{reaction_side_pattern})$"

    while True:
        reaction = input("Enter reaction formula:")
        reaction = reaction.rstrip()
        print(reaction, "\t", reaction_pattern)

        try:
            found = re.search(reaction_pattern, reaction)
            print(0, "\t found ", found)
            LH, RH = found.group("LH"), found.group("RH")
            print(1, "\t LH: ", LH, " and RH: ", RH)
            res = []
            for side in [LH, RH]:
                molecules = extract_molecules(molecule_pattern, side)
                molecules_data = []
                for molecule in molecules:
                   counts = extract_atom_data(molecule)
                   molecules_data.append(counts)
                res.append(molecules_data)
            
            return res
        except:
            print("Could not parse expression try again!")


if __name__ == '__main__':
    res = take_input()
    print(res)