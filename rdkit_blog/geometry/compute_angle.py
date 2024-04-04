from rdkit import Chem
from rdkit.Chem import rdMolTransforms

def get_bond_type(bond):
    if bond is None:
        return None
    elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
        return '-'
    elif bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
        return '='
    elif bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
        return '#'
    else:
        return '?'

def generate_angle_list(mol):
    angle_list = []
    conf = mol.GetConformer()
    for idx1 in range(mol.GetNumAtoms()):
        for idx2 in range(mol.GetNumAtoms()):
            for idx3 in range(mol.GetNumAtoms()):
                if idx1 < idx2 and idx2 < idx3:
                    angle = rdMolTransforms.GetAngleDeg(conf, idx1, idx2, idx3)
                    atom1, atom2, atom3 = mol.GetAtomWithIdx(idx1), mol.GetAtomWithIdx(idx2), mol.GetAtomWithIdx(idx3)
                    bond1_type = get_bond_type(mol.GetBondBetweenAtoms(idx1, idx2))
                    bond2_type = get_bond_type(mol.GetBondBetweenAtoms(idx2, idx3))
                    if bond1_type and bond2_type:
                        angle_type = f"{atom1.GetSymbol()}{bond1_type}{atom2.GetSymbol()}{bond2_type}{atom3.GetSymbol()}"
                        angle_list.append((angle_type, (idx1, idx2, idx3), angle))
    return angle_list

if __name__ == '__main__':

    mol = Chem.MolFromMolFile('./druglike_property/test.sdf')
    Chem.Kekulize(mol)

    angle_list = generate_angle_list(mol)
    print(angle_list)