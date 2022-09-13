from rdkit import Chem 
from rdkit.Chem import AllChem

def get_pocket(lig_mol, rec_mol, max_dist=8):
    lig_coords = lig_mol.GetConformer().GetPositions()
    rec_coords = rec_mol.GetConformer().GetPositions()
    dist = sp.spatial.distance.cdist(lig_coords, rec_coords)

    # indexes of atoms in rec_mol that are
    #   within max_dist of an atom in lig_mol
    pocket_atom_idxs = set(np.nonzero((dist < max_dist))[1])

    # determine pocket residues
    pocket_res_ids = set()
    for i in pocket_atom_idxs:
        atom = rec_mol.GetAtomWithIdx(int(i))
        res_id = get_rd_atom_res_id(atom)
        pocket_res_ids.add(res_id)

    # copy mol and delete atoms
    pkt_mol = rec_mol
    pkt_mol = Chem.RWMol(pkt_mol)
    for atom in list(pkt_mol.GetAtoms()):
        res_id = get_rd_atom_res_id(atom)
        if res_id not in pocket_res_ids:
            pkt_mol.RemoveAtom(atom.GetIdx())

    Chem.SanitizeMol(pkt_mol)
    return pkt_mol

rd_mol = Chem.RWMol(rd_mol)
rec_mol = Chem.MolFromPDBFile(os.path.join(data_root, mol_dicts['rec_src']), sanitize=True)
rec_mol = get_pocket(rd_mol, rec_mol)
uff_mol = Chem.CombineMols(rec_mol, rd_mol)
try:
    Chem.SanitizeMol(uff_mol)
except Chem.AtomValenceException:
    print('Invalid valence')
except (Chem.AtomKekulizeException, Chem.KekulizeException):
    print('Failed to kekulize')
try:
    uff = AllChem.UFFGetMoleculeForceField(
                    uff_mol, confId=0, ignoreInterfragInteractions=False
                )
    uff.Initialize()
    for i in range(rec_mol.GetNumAtoms()): # Fix the rec atoms
        uff.AddFixedPoint(i)
    converged = False
    n_iters=200
    n_tries=2
    while n_tries > 0 and not converged:
        print('.', end='', flush=True)
        converged = not uff.Minimize(maxIts=n_iters)
        n_tries -= 1
    print(flush=True)
    print("Performed UFF with binding site...")
except:
    print('Skip UFF...')
coords = uff_mol.GetConformer().GetPositions()
rd_conf = rd_mol.GetConformer()
for i, xyz in enumerate(coords[-rd_mol.GetNumAtoms():]):
    rd_conf.SetAtomPosition(i, xyz)

#得到的rd_conf 就是优化好的构象
