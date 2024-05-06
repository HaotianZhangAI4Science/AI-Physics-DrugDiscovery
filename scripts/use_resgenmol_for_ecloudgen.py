# /home/haotian/molecules_confs/Protein_test/Res2mol/val/docking_analysis.ipynb
# transfer the top 10 molecules from resgen_epo57 to resgen_epo57_top10
import os
resgen_old_base = './resgen_epo57'
resgen_new_base = './resgen_epo57_top10'
target_keys = list(mols.keys())
for target_name in target_keys:
    try:
        if len(mols[target_name]['gen_docking']) > 10:
            old_target = osp.join(resgen_old_base,target_name)
            new_target = osp.join(resgen_new_base,target_name)
            os.makedirs(new_target, exist_ok=True)
            ori_sdf_file = max(glob(osp.join(old_target, '*.sdf')), key=len)
            pdb_file = min(glob(osp.join(old_target, '*.pdb')), key=len)
            shutil.copy(pdb_file, new_target)
            shutil.copy(ori_sdf_file, new_target)
            write_sdf(mols[target_name]['gen_docking'][:10], osp.join(new_target, 'resgen_top10_mols.sdf'))
    except Exception as e:
        print(e)

# /home/haotian/Molecule_Generation/MG/ECloudGen_ELIP/generate_lig_ecloud.ipynb
# generate ligand ecloud for resgen_top10_mols.sdf
targets = glob('/home/haotian/Molecule_Generation/MG/backupECloud/EcloudGen-COATI/generation/resgen_epo57_top10/*')
for target in tqdm(targets):
    sdf_file = glob(f'{target}/resgen_top10_mols.sdf')[0]
    gen_mols = read_sdf(sdf_file)
    lig_density_list = []
    saved = f'{target}/resgen_top10.npy'
    if os.path.exists(saved):
        continue
    for mol in gen_mols:
        try:
            lig_density = get_ligecloud(mol,calculater, protocol(32))
            lig_density_list.append(lig_density)
        except Exception as e:
            print(e)
    np.save(f'{target}/resgen_top10.npy',lig_density_list)