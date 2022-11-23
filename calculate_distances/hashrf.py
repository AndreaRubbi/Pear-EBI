import subprocess
import numpy as np
import pandas as pd
import os

# Set the value Display variable
os.environ.setdefault('DISPLAY', ':0.0')

# Simple function to perform bash operations 
def bash_command(cmd): 
    subprocess.run(['/bin/bash', '-c', cmd])
    return 0

# ──────────────────────────────── RUNNING HASHRF RETURNING CLEANED OUTPUT ─────
def hashrf(file, n_trees, output_file):
    try: path = (os.path.dirname(__file__))
    except: path = '.'
    cmd = '{path}/HashRF/hashrf {file} {n_trees} | tail -{tail} | head -{n_trees} > {output}'.format(path = path,
        file= file, tail = n_trees + 3, n_trees= n_trees, output= output_file, verbose = 1)
        
    # Runs command
    bash_command(cmd) 
    
    # Reads output file and creates numpy array
    # which is rused to create a pandas dataframe
    with open(output_file, 'r') as out:
            distance_matrix = pd.DataFrame(np.loadtxt(out))
            out.close()
    
    distance_matrix.to_csv(output_file)
    return distance_matrix



        
        
        
        
        
        
        
        
    
        










