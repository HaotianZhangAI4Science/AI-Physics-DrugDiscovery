from rdkit import Chem
from rdkit.Chem import AllChem
import copy
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import AllChem, TorsionFingerprints


def GetEnergy(Mol, ConfID = None, FF_type='MMFF'):
    '''
    Calculate the molecular energy by empirical force field
    '''
    Status = True
    Energy = None

    if ConfID is None:
        ConfID = -1
    
    if FF_type == "UFF":
        UFFMoleculeForcefield = AllChem.UFFGetMoleculeForceField(Mol, confId = ConfID)
        if UFFMoleculeForcefield is None:
            Status = False
        else:
            Energy = UFFMoleculeForcefield.CalcEnergy()
    elif FF_type == "MMFF":
        #In the RDKit, MMFF has two options, “MMFF94” / “MMFF94s”.
        MMFFMoleculeProperties = AllChem.MMFFGetMoleculeProperties(Mol)
        MMFFMoleculeForcefield = AllChem.MMFFGetMoleculeForceField(Mol, MMFFMoleculeProperties, confId = ConfID)
        if MMFFMoleculeForcefield is None:
            Status = False
        else:
            Energy = MMFFMoleculeForcefield.CalcEnergy()
    else:
        print('Only support the UFF or MMFF!')
    
    return (Status, Energy)


def rigid_torsion_scan(mol, start_angle, end_angle, step_angle, torsion_ids, addHs=True, \
    gen_init_conf=False, FF_option='MMFF'):
    '''
    Perform the rigid torsion scan
    torsion_id: a list which contains 4 int to specify the dihedral angle
    '''
    idx_1, idx_2, idx_3, idx_4 = torsion_ids
    dihedral_energies = []
    angles = [Angle for Angle in range(start_angle, end_angle, step_angle)]
    if addHs:
        mol = Chem.AddHs(mol)
    if gen_init_conf:
        AllChem.EmbedMolecule(mol)
        if FF_option == 'MMFF':
            AllChem.MMFFOptimizeMolecule(mol)
        elif FF_option == 'UFF':
            AllChem.UFFOptimizeMolecule(mol)
        else:
            print('No optimization since the FF option is not supported')
    mol_ = copy.deepcopy(mol)
    for dihedral_degree in angles:
        Torsionmol = Chem.Mol(mol_)
        conformer = Torsionmol.GetConformer(0)
        rdMolTransforms.SetDihedralDeg(conformer, idx_1, idx_2, idx_3, idx_4, dihedral_degree)
        CalcStatus, Energy = GetEnergy(Torsionmol,FF_type=FF_option)
        dihedral_energies.append(Energy)
    return angles, dihedral_energies


def print_torsions(mol):
    nonring, ring = TorsionFingerprints.CalculateTorsionLists(mol)
    conf = mol.GetConformer(id=0)
    tups = [atoms[0] for atoms, ang in nonring]
    degs = [Chem.rdMolTransforms.GetDihedralDeg(conf, *tup) for tup in tups]
    print(degs)