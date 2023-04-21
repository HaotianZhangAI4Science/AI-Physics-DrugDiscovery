```python
conda install -c conda-forge pymol-open-source
pip install freesasa # installation from conda can not work for me. 
mamba install -c conda-forge fpocket # fpocket has no api, it can just be called by terminal. e.g. fpocket -f af_pocket.pdb
conda install mmseq2 -c bioconda
```

mmseqs compute sequence similarity example:

```python
# firstly, create mmseqs database
mmseqs createdb train.fasta train_db
mmseqs createdb test.fasta test_db
# secondly, search seqs similarity
mmseqs search test_db train_db
# finally, convert the search results into a tabular format
mmseqs convertalis test_db train_db search_result search_result.tsv
# read search_result.tsv
train_set_ids = {}  # List of unique identifiers for the train set
test_set_ids = {}  # List of unique identifiers for the test set

with open("./tmp/search_result.tsv") as f:
    for line in f:
        query, target, similarity, _,_, _, _, _, _, _, _, _ = line.strip().split("\t")
        similarity = float(similarity)
        try:
            test_set_ids[query].append(similarity)
        except:
            test_set_ids[query] = []
            test_set_ids[query].append(similarity)
```

parKVFiner is another cavity detection method.

[LBC-LNBio/parKVFinder: parKVFinder: thread-level parallel KVFinder (github.com)](https://github.com/LBC-LNBio/parKVFinder)

## pdb_utils

#### pdb_to_fasta

#### read_fasta_file

#### pdb_to_fasta_parallel

#### align_and_rmsd

#### residues_saver

#### pocket_trunction

#### extract_alphafold_pocket

#### calculate_pocket_sasa

####  read_pocket_info

#### compute_convex_hull_volume

#### find_diSContact  (in disulfide.py)

return the sulfide-bonded atoms



## alphafold_utils

#### pdb2uniprot

#### download_alphafold_prediction
