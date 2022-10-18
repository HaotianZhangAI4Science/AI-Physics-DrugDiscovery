rdkit rotate molecules https://iwatobipen.wordpress.com/2019/02/27/rotate-molecule-and-visualize-it-rdkit-chemoinformatics/
rdkit rtz https://github.com/wutobias/r2z
https://chemistry.stackexchange.com/questions/136930/computational-chemistry-software-to-generate-files-of-cartesian-coordinates-and
rdkit切分fragment
https://zhuanlan.zhihu.com/p/390678258
rdkit显示键的编号
https://zhuanlan.zhihu.com/p/163827512
二面角扫描
https://pschmidtke.github.io/blog/torsion/dihedral/oss/opensource/rdkit/xtb/energy/2021/02/16/torsion-angle-scans-xtb.html
https://github.com/vermasrijan/redial-2020/blob/5cfd09d5645301654373cb9008ded1b2552bd920/mayachemtools/bin/RDKitPerformTorsionScan.py#L588
fragment_utils
https://github.com/deepchem/deepchem/blob/9664aeab16940e1ee46030cac08d90e12aa24751/deepchem/utils/test/test_fragment_utils.py

不同rdkit版本的容错率还不一样，如果直接从mol文件当中读取分子，比如
mol = Chem.MolFromMolFile(lig_path) 的话，很多时候会报valence的错误
这个时候mol = Chem.MolFromMolFile(lig_path,sanitize=False) 就可以继续操作了

读取sdf文件
suppl = Chem.SDMolSupplier(lig_path, sanitize=False)
mols = [mol for mol in suppl if mol]	
mol = mols[0]

读取mol2文件
mol = Chem.MolFromMol2File(mol2_path, sanitize=False)

屏蔽Rdkit的Error输出
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
smi = 'CO(C)C'
mol = Chem.MolFromSmiles(smi)
try:
    Chem.SanitizeMol(mol)
except:
    print('OK')
   
# rdkit内部会有拓扑匹配的事情，所以要非常小心地处理拓扑问题。拿到新的分子最好的方式就是先去H然后再加氢
def cannaical_mol(mol, 3d=True, addhs=True):
    mol_ = copy.deepcopy(mol)
    mol_ = Chem.RemoveAllHs(mol_)
    smi = Chem.MolToSmiles(mol_)
    if 3D is True and addhs is True:
        mol_ = Chem.AddHs(mol_, addCoords=True)
    if 3D is not True and addhs is True:
        mol_ = Chem.AddHs(mol_)
    return mol_, smi
    return mol_, smi

def rm_radical(mol):
    mol = copy.deepcopy(mol)
    for atom in mol.GetAtoms():
        atom.SetNumRadicalElectrons(0)
    return mol
