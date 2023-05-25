from Bio.PDB import PDBParser, NeighborSearch, PDBIO, Select, PPBuilder
from pymol import cmd
import os.path as osp
import os
from rdkit import Chem
# conda install -c conda-forge pymol-open-source
from Bio import PDB
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

def pdb_to_fasta(pdb_filename, output_fasta):
    pp_builder = PPBuilder()
    parser = PDBParser(QUIET=True)
    records = []


    structure = parser.get_structure(os.path.splitext(os.path.basename(pdb_filename))[0], pdb_filename)
        
    all_chains_sequence = ""
    for model in structure:
        for chain in model:
            pp = pp_builder.build_peptides(chain)
            sequence = "".join(str(p.get_sequence()) for p in pp)
            all_chains_sequence += sequence
                
    seq_record = SeqRecord(Seq(all_chains_sequence), id=structure.id, description=f"All chains from {structure.id}")

    SeqIO.write([seq_record], output_fasta, "fasta")


def read_fasta_file(file_path):
    with open(file_path, "r") as file:
        lines = file.readlines()
    sequences = {}
    current_key = ""
    for line in lines:
        if line.startswith(">"):
            current_key = line.strip()
            sequences[current_key] = ""
        else:
            sequences[current_key] += line.strip()
    return sequences


def pdb_to_fasta_parallel(pdb_filenames, output_fasta, num_workers=4): 
    '''
    e.g.: pdb_to_fasta_parallel(train_pdbfiles, "train.fasta", num_workers=12)

    '''
    def process_pdb_file(pdb_filename):
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure(os.path.basename(pdb_filename), pdb_filename)
        name = os.path.basename(pdb_filename)
        sequences = []
        seq = ""
        for chain in structure.get_chains():
            for residue in chain.get_residues():
                res_id = residue.get_resname()
                try:
                    seq += PDB.Polypeptide.three_to_one(res_id)
                except KeyError:
                    continue

        sequences.append(f">{name}\n{seq}")

        return sequences

    all_sequences = []

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(process_pdb_file, pdb_filename): pdb_filename for pdb_filename in pdb_filenames}

        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing PDB files"):
            sequences = future.result()
            all_sequences.extend(sequences)

    with open(output_fasta, "w") as output_file:
        output_file.write("\n".join(all_sequences))




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
    pocket_residues = {res for coord in ligand_coords for res in ns.search(coord, threshold)}
    
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

    return pocket_residues, pocket_residues_alphafold
    
def remove_hetatom(pdb_file, out_file=None):

    protein_dir = os.path.dirname(pdb_file)
    pdb_name = os.path.basename(pdb_file).split('.')[0]

    temp_path_0 = os.path.join(protein_dir, f'{pdb_name}_temp_0.pdb')
    temp_path_1 = os.path.join(protein_dir, f'{pdb_name}_temp_1.pdb')

    if out_file is None:
        out_protein_path = os.path.join(protein_dir, f'{pdb_name}_out.pdb')
    else:
        out_protein_path = out_file

    subprocess.run(f'pdb_selmodel -1 {pdb_file} > {temp_path_0}', shell=True)
    subprocess.run(f'pdb_delelem -H {temp_path_0} > {temp_path_1}', shell=True)
    subprocess.run(f'pdb_delhetatm {temp_path_1} > {out_protein_path}', shell=True)

    os.remove(temp_path_0)
    os.remove(temp_path_1)


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
