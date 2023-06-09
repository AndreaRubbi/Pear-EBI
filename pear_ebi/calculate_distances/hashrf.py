import os
import subprocess
import sys

import numpy as np
import pandas as pd

# Set the value Display variable
os.environ.setdefault("DISPLAY", ":0.0")

# Simple function to perform bash operations


def bash_command(cmd):
    """Executes bash command in subprocess

    Args:
        cmd (str): bash command to be runned in subprocess

    Returns:
        0: returns 0 if everything's alright
    """
    if os.name == "nt":
        sys.exit("PEAR isn't compatble with Windows yet")
        try:
            subprocess.run(["C:\\cygwin64\\bin\\bash.exe", "-lc", cmd])
        except:
            sys.exit("Windows platforms need cygwin64/cygwin to run this subprocess")
    elif os.name == "posix":
        subprocess.run(
            ["/bin/bash", "-c", cmd], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL
        )
    else:
        sys.exit("Could not verify OS")
    return 0


# ──────────────────────────────── RUNNING HASHRF RETURNING CLEANED OUTPUT ─────
# ? From: "A RANDOMIZED ALGORITHM FOR COMPARING SETS OF PHYLOGENETIC TREES"
# ? by Seung-Jin Sul & Tiffani L. Williams
# ? https://code.google.com/archive/p/hashrf/
# ──────────────────────────────────────────────────────────────────────────────


def hashrf(file, n_trees, output_file):
    """Computes unweighted Robison Foulds distances

    Args:
        file (str): name of input file with phylogenetic trees in newick format
        n_trees (int): number of trees in file
        output_file (str): name of output file that will contain the distance matrix

    Returns:
        distance_matrix (pandas.DataFrame): computed distance matrix
    """
    try:
        path = os.path.dirname(os.path.abspath(__file__))
    except:
        path = "."
    cmd = "{path}hashrf {file} {n_trees} -p matrix -o {output} 1> /dev/null 2> ./hashrf_err.txt".format(
        path=os.path.join(path, "HashRF", ""),
        file=file,
        n_trees=n_trees,
        output=output_file,
        verbose=1,
    ).replace(
        "\\", "/"
    )

    # Runs command
    bash_command(cmd)

    try:
        with open("hashrf_err.txt", "r") as err_mess:
            err = err_mess.read()
            err_mess.close()
            bash_command("rm hashrf_err.txt")
            if len(err) > 0:
                sys.exit(err)
    except:
        pass

    # Reads output file and creates numpy array
    # which is used to create a pandas dataframe
    try:
        with open(output_file, "r") as out:
            distance_matrix = pd.read_csv(out, index_col=None, header=None, sep=" ")
            out.close()

        distance_matrix.drop(distance_matrix.columns[-1], axis=1, inplace=True)
        distance_matrix.to_csv(output_file, header=False, index=False)
        return distance_matrix.values
    except:
        sys.exit("hashrf failed!")


# HashRF calculating weighted RF distances
def hashrf_weighted(file, n_trees, output_file):
    """Computes weighted Robison Foulds distances

    Args:
        file (str): name of input file with phylogenetic trees in newick format
        n_trees (int): number of trees in file
        output_file (str): name of output file that will contain the distance matrix

    Returns:
        distance_matrix (pandas.DataFrame): computed distance matrix
    """
    try:
        path = os.path.dirname(__file__)
    except:
        path = "."
    cmd = "{path}/HashRF/hashrf {file} {n_trees} -p matrix -o {output} -w 1> /dev/null 2> ./hashrf_err.txt".format(
        path=path, file=file, n_trees=n_trees, output=output_file, verbose=1
    )

    # Runs command
    bash_command(cmd)

    try:
        with open("hashrf_err.txt", "r") as err_mess:
            err = err_mess.read()
            err_mess.close()
            bash_command("rm hashrf_err.txt")
            if len(err) > 0:
                sys.exit(err)
    except:
        pass

    bash_command(
        f"tr -s \" \" < {output_file} | sed 's/^[ \t]*//' > ./tmp_file && cat ./tmp_file > {output_file} && rm ./tmp_file"
    )  # removes double spaces in file and spaces at the beginning of lines

    # Reads output file and creates numpy array
    # which is rused to create a pandas dataframe
    try:
        with open(output_file, "r") as out:
            distance_matrix = pd.read_csv(out, index_col=None, header=None, sep=" ")
            out.close()

        distance_matrix.drop(distance_matrix.columns[-1], axis=1, inplace=True)
        distance_matrix.to_csv(output_file, header=False, index=False)
        return distance_matrix.values
    except:
        sys.exit("hashrf failed!")
