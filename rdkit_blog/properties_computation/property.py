from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
import os.path as osp
import mdtraj as md
import numpy as np
import os


def hbd(mol):
    return rdMolDescriptors.CalcNumHBA(mol)
def hbd(mol):
    return rdMolDescriptors.CalcNumHBD(mol)
def tpsa(mol):
    return rdMolDescriptors.CalcTPSA(mol)

import mdtraj as md
import numpy as np

def sasa(mol):
    pdb_file = './tmp'
    Chem.MolToPDBFile(mol, pdb_file)
    sasa = compute_sasa(pdb_file)
    os.remove(pdb_file)
    return sasa

def compute_sasa(pdb_file):
    # Load the PDB file
    traj = md.load(pdb_file)

    # Compute the SASA
    sasa = md.shrake_rupley(traj, mode='residue')

    # Get the total SASA
    total_sasa = np.sum(sasa)

    return total_sasa