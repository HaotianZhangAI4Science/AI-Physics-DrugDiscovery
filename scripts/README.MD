### Python Command 

We can automatically run the scripts in the python file instead in the bash file. 

```python
import subprocess 

command = f'python -u delete.py --surf_path {ply_file} --frag_path {frag_file} --check_point {ckpt} --outdir ./{out_dir} --suboutdir {suboutdir} '
print(command)
result = subprocess.run(command, shell=True, capture_output=True, text=True)

if result.returncode == 0:
    print('executed successfully.')
    print('Output:')
    print(result.stdout)
    print('consumed time: ',time.time()-start_time)
else:
    print('execution failed.')
    print('Error:')
    print(result.stderr)
```



### Python multiprocess execution

We use the process of pdbfixer as an example 

```python
from multiprocessing import Pool
import subprocess
from glob import glob
import os.path as osp
from tqdm import tqdm

# define a function that executes your command 
def execute_command(pdb_file):
    pdb_out_file = pdb_file[:-4] + '_added.pdb'
    command = f'pdbfixer {pdb_file} --output {pdb_out_file} --replace-nonstandard --add-atoms=all --add-residues'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        return f'Execution failed for {pdb_file}. Error: {result.stderr}'
    return f'Execution successful for {pdb_file}'

# Use glob to find all directories in 'crossdock' and then find the first .pdb file in each
targets = glob('./crossdock/*')
pdb_files = []
for target in targets:
    pdb_files.extend(glob(osp.join(target, '*.pdb')))
    
# Create a pool of worker processes, with the number of processes set to the desired level
# This is often set to the number of CPU cores on the machine
pool = Pool(processes=24)

# Use tqdm to create a progress bar for the execution of the pool
results = list(tqdm(pool.imap(execute_command, pdb_files), total=len(pdb_files)))

# Close the pool to prevent any more tasks from being submitted to the pool
pool.close()

# Wait for the worker processes to finish
pool.join()

# Print out the results of the command executions
for result in results:
    print(result)

```

