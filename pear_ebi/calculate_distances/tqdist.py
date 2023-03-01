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
        try:
            subprocess.run(
                ["C:\\cygwin64\\bin\\bash.exe", "-l", cmd],
                stderr=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL,
            )
        except:
            try:
                subprocess.run(
                    ["C:\\cygwin\\bin\\bash.exe", "-l", cmd],
                    stderr=subprocess.DEVNULL,
                    stdout=subprocess.DEVNULL,
                )
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
# ? From: "tqDist: a library for computing the quartet and triplet distances between binary or general trees"
# ? by A. Sand, C. N. S. Pedersen et al - 2014
# ? https://www.birc.au.dk/~cstorm/software/tqdist/
# ──────────────────────────────────────────────────────────────────────────────


def quartet(file, n_trees, output_file):
    """Computes quartet distances

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
    cmd = "{path}/tqDist/bin/all_pairs_quartet_dist {file} {output} 1> /dev/null 2> ./tqdist_err.txt".format(
        path=path, file=file, output=output_file, verbose=1
    )

    # Runs command
    bash_command(cmd)
    bash_command(
        f"tr -s \" \" < {output_file} | sed 's/\s/_/g' > ./tmp_file && cat ./tmp_file > {output_file} && rm ./tmp_file"
    )

    # Reads output file and creates numpy array
    # which is used to create a pandas dataframe
    try:
        with open(output_file, "r") as out:
            distance_matrix = np.zeros((n_trees, n_trees))
            for i, line in enumerate(out.readlines()):
                values = line.strip().split("_")[:-1]
                distance_matrix[i, : i + 1] = values
            out.close()

        distance_matrix_upper = distance_matrix.transpose()
        distance_matrix = pd.DataFrame(distance_matrix + distance_matrix_upper)
        distance_matrix.to_csv(output_file)
        return distance_matrix
    except:
        print("quartet failed! check errors in hashrf_err.txt")


def triplet(file, n_trees, output_file):
    """Computes triplet distances

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
    cmd = "{path}/tqDist/bin/all_pairs_quartet_dist {file} {output} 1> /dev/null 2> ./tqdist_err.txt".format(
        path=path, file=file, output=output_file, verbose=1
    )

    # Runs command
    bash_command(cmd)
    bash_command(
        f"tr -s \" \" < {output_file} | sed 's/\s/_/g' > ./tmp_file && cat ./tmp_file > {output_file} && rm ./tmp_file"
    )

    # Reads output file and creates numpy array
    # which is rused to create a pandas dataframe
    try:
        with open(output_file, "r") as out:
            distance_matrix = np.zeros((n_trees, n_trees))
            for i, line in enumerate(out.readlines()):
                values = line.strip().split("_")[:-1]
                distance_matrix[i, : i + 1] = values
            out.close()

        distance_matrix_upper = distance_matrix.transpose()
        distance_matrix = pd.DataFrame(distance_matrix + distance_matrix_upper)
        distance_matrix.to_csv(output_file)
        return distance_matrix
    except:
        print("triplet failed! check errors in hashrf_err.txt")
