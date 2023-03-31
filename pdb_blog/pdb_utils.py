from Bio.PDB import PDBParser, NeighborSearch, PDBIO, Select
from pymol import cmd
import os.path as osp
from rdkit import Chem
# conda install -c conda-forge pymol-open-source

def align_and_rmsd(src_file, tgt_file, saved_file=None):
    '''
    align the src_file to the tgt_file, and saved src_file to the saved_file

    '''
    cmd.load(src_file, "src_structure")
    cmd.load(tgt_file, "tgt_structure")

    # Align the structures
    align_results = cmd.align("src_structure", "tgt_structure")

    # The 1st element of align_results is the RMSD value
    rmsd = align_results[0]
    # Save the aligned structures
    if saved_file is not None:
        cmd.save(saved_file, "src_structure")
        print('PyMol saved: ', saved_file)

    # Clean up
    cmd.delete("src_structure")
    cmd.delete("tgt_structure")

    return rmsd



class PocketResidueSelect(Select):
    def __init__(self, residues):
        self.residues = residues

    def accept_residue(self, residue):
        return residue in self.residues


def residues_saver(structure,residues, out_name, verbose=0):
    pdbio = PDBIO()
    pdbio.set_structure(structure)
    pdbio.save(out_name, PocketResidueSelect(residues))
    if verbose:
        print('residues saved at {}'.format(out_name))

def load_ligand(ligand):
    
    ligand_name = osp.basename(ligand)
    suffix = ligand_name.split('.')[-1]
    if suffix == 'mol2':
        ligand_mol = Chem.MolFromMol2File(ligand)
    elif suffix == 'mol' or suffix == 'sdf':
        ligand_mol = Chem.MolFromMolFile(ligand)
    else:
        ligand_mol = ligand
    
    return ligand_mol


def pocket_trunction(pdb, ligand, threshold=10.0, save_name=None):
    '''
    trunct pocket centered at `ligand within `threshold
    args:
    pdb:       `PDB file or Biopython structure
    ligand:    `ligand file or rdkit mol
    threshold: `pocket radius
    '''
    # Load Ligand Structures 
    ligand_mol = load_ligand(ligand)
    
    ligand_conf = ligand_mol.GetConformer()
    ligand_coords = ligand_conf.GetPositions()

    # Load Protein Sructures
    if type(pdb) == str:
        parser = PDBParser()
        structure = parser.get_structure("crystal", pdb)
    else:
        structure = pdb
    # Extract pocket residues
    ns = NeighborSearch(list(structure.get_atoms()))
    pocket_residues = {res.parent for coord in ligand_coords for res in ns.search(coord, threshold)}
    
    # (Optional): save the pocket resides
    if save_name:
        residues_saver(structure, pocket_residues, save_name)
    
    return pocket_residues


def extract_alphafold_pocket(pocket_residues, alphafold_structure):

    # Get residue names and ids for the pocket residues
    pocket_residue_info = {(res.get_resname(), res.id[1]) for res in pocket_residues}

    # Extract pocket residues in the AlphaFold structure
    pocket_residues_alphafold = set()
    for chain in alphafold_structure.get_chains():
        for res in chain:
            if (res.get_resname(), res.id[1]) in pocket_residue_info:
                pocket_residues_alphafold.add(res)

    return pocket_residues_alphafold
    

def calculate_pocket_sasa(pocket_pdb_file):
    import freesasa
    # Load the structure from the PDB file
    structure = freesasa.Structure(pocket_pdb_file)

    # Compute the SASA
    result = freesasa.calc(structure)

    # Get the total SASA
    sasa = result.totalArea()

    return sasa

def read_pocket_info(pdb_file):
    pocket_coords = []
    monte_carlo_volume = None
    convex_hull_volume = None

    with open(pdb_file, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                pocket_coords.append([x, y, z])
            if "Pocket volume (Monte Carlo)" in line:
                monte_carlo_volume = float(line.split(":")[-1].strip())
            if "Pocket volume (convex hull)" in line:
                convex_hull_volume = float(line.split(":")[-1].strip())
    
    return np.mean(pocket_coords, axis=0), monte_carlo_volume, convex_hull_volume

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull

def compute_convex_hull_volume(coords):
    tri = Delaunay(coords)
    tetrahedra = coords[tri.simplices]

    volume = 0.0
    for tetrahedron in tetrahedra:
        volume += np.abs(np.linalg.det(tetrahedron[:-1] - tetrahedron[-1])) / 6

    return volume

def plot_convex_hull(coords1, coords2, volume1, volume2):
    hull1 = ConvexHull(coords1)
    hull2 = ConvexHull(coords2)

    fig = plt.figure(figsize=(12, 6))
    ax1 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122, projection='3d')

    for ax, hull, volume, title in zip([ax1, ax2], [hull1, hull2], [volume1, volume2], ['Pocket 1', 'Pocket 2']):
        for simplex in hull.simplices:
            poly = Poly3DCollection([coords1[simplex]], alpha=0.5)
            poly.set_facecolor('blue')
            ax.add_collection3d(poly)

        ax.set_title(f'{title} - Volume: {volume:.2f}')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_xlim([-50, 50])
        ax.set_ylim([-50, 50])
        ax.set_zlim([-50, 50])
        ax.view_init(elev=20, azim=-35)

    plt.show()


def get_pdb_coords(pdbfile):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('pocket', pdbfile)
    coords = np.array([atom.get_coord() for atom in structure.get_atoms()])
    return coords

# compute and plot convex hull
# coords1 = np.random.rand(30, 3) * 100 - 50  # Generate random coordinates for pocket 1
# coords2 = np.random.rand(30, 3) * 100 - 50  # Generate random coordinates for pocket 2
# coords1 = get_pdb_coords(crystal)
# coords2 = get_pdb_coords(af)
# volume1 = compute_convex_hull_volume(coords1)  # Replace with your function to compute the volume
# volume2 = compute_convex_hull_volume(coords2)  # Replace with your function to compute the volume

# plot_convex_hull(coords1, coords2, volume1, volume2)

# Align and map 
# from Bio import pairwise2
# from Bio.SubsMat import MatrixInfo

# def extract_pocket_residues(structure, alphafold_structure, ligand_coords, distance_threshold):
#     ns = NeighborSearch(list(structure.get_atoms()))

#     # Find pocket residues in the crystal structure
#     pocket_residues = {res for coord in ligand_coords for res in ns.search(coord, distance_threshold)}

#     # Align sequences
#     crystal_sequence = "".join([res.get_resname() for res in structure.get_residues()])
#     alphafold_sequence = "".join([res.get_resname() for res in alphafold_structure.get_residues()])
#     alignment = pairwise2.align.globalds(crystal_sequence, alphafold_sequence, MatrixInfo.blosum62, -10, -0.5)[0]

#     # Get residue indices of pocket residues in the crystal structure
#     pocket_residue_indices_crystal = [i for i, res in enumerate(structure.get_residues()) if res in pocket_residues]

#     # Map residue indices from the crystal structure to the AlphaFold structure using the sequence alignment
#     crystal_to_alphafold_map = {}
#     crystal_index, alphafold_index = 0, 0
#     for (i, j) in zip(alignment[0], alignment[1]):
#         if i != "-":
#             crystal_index += 1
#         if j != "-":
#             alphafold_index += 1
#         if i != "-" and j != "-":
#             crystal_to_alphafold_map[crystal_index - 1] = alphafold_index - 1

#     pocket_residue_indices_alphafold = [crystal_to_alphafold_map[i] for i in pocket_residue_indices_crystal]

#     # Extract pocket residues in the AlphaFold structure
#     pocket_residues_alphafold = [res for i, res in enumerate(alphafold_structure.get_residues()) if i in pocket_residue_indices_alphafold]

#     return pocket_residues, pocket_residues_alphafold
