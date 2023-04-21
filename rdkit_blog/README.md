# RDKit Toolkit and tutorial

### highlight_scaffold

Draw the picture like the following:

<div align=center>
<img src='./highlight_scaffold/output.png'width="50%"height="50%"align=center />
</div>
### properties_computation

#### hba(mol)

#### hdb(mol)

#### tpsa(mol)

#### sasa(mol)

#### compute_sasa(pdb_file)



```python
def mol_with_atom_index(mol):
    atoms = mol.GetNumAtoms()
    for idx in range(atoms):
        mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(mol.GetAtomWithIdx(idx).GetIdx()))
    return mol

except Exception as e:
    print(f'error during 3D generation -- {e}')
    continue
```

