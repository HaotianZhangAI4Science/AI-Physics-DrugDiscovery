#download pdb
import subprocess
import shutil
for i in range(len(pdb_id)):
    command = 'wget http://www.rcsb.org/pdb/files/{}.pdb'.format(pdb_id[i])
        
    proc = subprocess.Popen(
            command, 
            shell=True, 
            stdin=subprocess.PIPE, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE
        )
    proc.communicate()
    shutil.move('{}.pdb'.format(pdb_id[i]), './pdb')
