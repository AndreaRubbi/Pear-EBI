import subprocess
import numpy as np
import pandas as pd
import os

# Set the value Display variable
os.environ.setdefault('DISPLAY', ':0.0')

# Simple function to perform bash operations 
def bash_command(cmd): 
    subprocess.run(['/bin/bash', '-c', cmd], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    return 0

# ──────────────────────────────── RUNNING HASHRF RETURNING CLEANED OUTPUT ─────
def hashrf(file, n_trees, output_file):
    try: path = (os.path.dirname(__file__))
    except: path = '.'
    cmd = '{path}/HashRF/hashrf {file} {n_trees} -p matrix -o {output} 1> /dev/null 2> ./hashrf_err.txt'.format(path = path,
        file= file, n_trees= n_trees, output= output_file, verbose = 1)
    
    # Runs command
    bash_command(cmd) 
    
    # Reads output file and creates numpy array
    # which is rused to create a pandas dataframe
    try:
        with open(output_file, 'r') as out:
                distance_matrix = pd.read_csv(out, index_col=None, header=None, sep=' ')
                out.close()
        
        distance_matrix.drop(distance_matrix.columns[-1], axis=1, inplace=True)
        distance_matrix.to_csv(output_file)
        return distance_matrix
    except: print("hashrf failed! check errors in hashrf_err.txt")



        
        
        
        
        
        
        
        
    
        










