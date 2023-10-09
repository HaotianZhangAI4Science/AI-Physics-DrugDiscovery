# RDKit Toolkit and tutorial

### highlight_scaffold

Draw the picture like the following:

<div align=center>
<img src='./highlight_scaffold/output.png'width="50%"height="50%"align=center />
</div>
func draw_align_mols in ./highlight_scaffold/draw_align.py

<div align=center>
<img src='./highlight_scaffold/output2.png'width="50%"height="50%"align=center />
</div>

### properties_computation

#### hba(mol)

#### hdb(mol)

#### tpsa(mol)

#### sasa(mol)

#### compute_sasa(pdb_file)

### Grid Representation

Do you want to make the representation as follows? Please find the function in the grid_repre repo. It contains grid representation computation with htmd package and visualization file style, i.e., cube file. 

<div align=center>
<img src='./grid_repre/grid.png'width="50%"height="50%"align=center />
</div>





```python
def mol_with_atom_index(mol):
    atoms = mol.GetNumAtoms()
    for idx in range(atoms):
        mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(mol.GetAtomWithIdx(idx).GetIdx()))
    return mol

from rdkit import Chem

def mol2svg(mol, file_name):
    mc = Chem.Mol(moltoBinary())
    if not mc.GetNumConformers():
        Chem.rdDepictor.Compute2DCoords(mc)
    drawer = Chem.Draw.rdMolDraw2D.MolDraw2DSVG(300,300)
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    svg_cleaned = svg.replace('svg:', '')
    with open(file_name, 'w') as f:
        f.write(str(svg_cleaned.data))
    

```

