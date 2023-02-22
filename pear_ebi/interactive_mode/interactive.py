import os
import sys

import rich
from rich import print

# getting the name of the directory
current = os.path.dirname(os.path.realpath(__file__))

# Getting the parent directory name
parent = os.path.dirname(current)

# adding the parent directory to
# the sys.path.
sys.path.append(parent)


def usage():
    print("[bold green]PEAR | Interactive Mode")
    print(
        "[green]⣿⣿⣿⣿⣿⣿⣿⣿⠿⣟⠉⡿⠿⣿⣿⣿⣿⣿⣿⣿[bold white] -- Controls --",
        "[green]⣿⣿⣿⣿⡿⣿⢉⢳⠴⣞⠉⡷⢥⡏⡙⡿⢿⣿⣿⣿[bold white] 1 --> see status",
        "[green]⣿⣿⡋⢻⡤⣼⠉⢯⡆⣞⠙⣧⣢⠏⠪⣣⢦⡛⠹⣿[bold white] 2 --> calculate distances",
        "[green]⣿⣿⠓⢻⣄⣼⠋⢷⡠⡽⠚⣉⣤⡞⢚⢦⢢⠟⠹⣿[bold white] 3 --> embed distances",
        "[green]⣿⣿⢓⢻⡄⡼⠗⢃⣂⡒⠻⣧⣂⡿⠚⣨⢨⡓⠻⣿[bold white] 4 --> plot embeddings",
        "[green]⣿⣿⠗⢎⢄⠂⠾⣯⣂⡽⢓⢆⢔⠐⠿⣇⢅⡗⠻⣿[bold white] 5 --> add set to collection",
        "[green]⣿⣿⠖⢯⡡⣹⠗⣤⣉⠛⠶⣏⢌⡿⠲⡌⢌⢞⠼⣿[bold white] 6 --> get subset",
        "[green]⣿⣿⣮⣾⡉⣹⠦⣞⡉⡽⠶⡌⠌⡮⠲⣏⢩⣷⣼⣿[bold white]",
        "[green]⣿⣿⣿⣿⣿⣿⣤⣞⢉⣳⠥⣏⠍⣧⣵⣿⣿⣿⣿⣿[bold white] 7 --> exit",
        "[green]⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣤⣿⣿⣿⣿⣿⣿⣿⣿⣿[bold white] 8 --> see list of controls",
        sep="\n",
    )

    return 0


def exit_program():
    print("[orange1]- Leaving PEAR -")
    exit()


def calculate_distances(SET):
    methods = ["hashrf", "hashrf_weighted", "days_RF", "quartet", "triplet"]
    while True:
        try:
            method = int(
                input(
                    "Method (1:hashrf - 2:weighted hashrf - 3:day's rf - 4:quartet - 5:triplet): "
                )
            )
        except ValueError:
            continue
        if method in (1, 2, 3, 4, 5):
            break
        else:
            rich.print("[bold orange1]Please select any available method")

    SET.calculate_distances(methods[method - 1])
    return 0


def embedding(SET):
    while True:
        try:
            method = int(input("Method (1:PCoA - 2:t-SNE): "))
        except ValueError:
            print("")
            continue
        if method in (1, 2):
            break
        else:
            print("[bold orange1]Please select either 1 or 2")

    while True:
        try:
            dimensions = int(input("Number of dimensions: "))
            break
        except ValueError:
            print("")
            continue

    if method == 1:
        SET.embed("pca", dimensions, quality=True)
    elif method == 2:
        SET.embed("tsne", dimensions, quality=True)
    return 0


def plotting(SET):
    methods = ["pca", "tsne"]
    while True:
        try:
            method = int(input("Method (1:PCoA - 2:t-SNE): "))
        except ValueError:
            print("")
            continue
        if method in (1, 2):
            break
        else:
            print("[bold orange1]Please select either 1 or 2")

    while True:
        try:
            dimensions = int(input("Number of dimensions: "))
        except ValueError:
            print("")
            continue
        if dimensions in (2, 3):
            break
        else:
            print("[bold orange1]Please select either 2 or 3")

    while True:
        show = input("Show figure (y/n): (default is no)  ")
        if show in ("y", "n", "Y", "N", "Yes", "No", "yes", "no", ""):
            break
        else:
            print("[bold orange1]Please select valid value")

    if dimensions == 2:
        fig = SET.plot_2D(methods[method - 1], save=True)
    elif dimensions == 3:
        fig = SET.plot_3D(methods[method - 1], save=True)

    if show in ("y", "Y", "Yes", "yes"):
        fig.show()
    return 0


def add_set():
    global SET
    while True:
        filename = input("Specify file name: ")
        try:
            open(filename, "r").close()
            break
        except:
            print("[bold orange1]Cannot find the specified file")
    return filename


def get_subset(SET):
    while True:
        while True:
            try:
                method = int(input("Method (1:syst/2:random/3:sequence): "))
                if method in (1, 2, 3):
                    break
            except ValueError:
                print("[bold orange1]Please select valid value")
        try:
            n_required = int(input("Number of elements in subset: "))
            if n_required >= SET.n_trees:
                print("Sample size larger or equal to population")
                return 1

            show = input("Show subset embedding (y/n): (default is no)  ")
            if show in ("y", "n", "Y", "N", "Yes", "No", "yes", "no", ""):
                break
        except ValueError:
            print("")
            continue

    fig1, fig2 = SET.get_subset(n_required, ["syst", "random", "sequence"][method - 1])
    if show in ("y", "Y", "Yes", "yes"):
        fig2.show()
    return 0


def interact(control):
    try:
        Actions = {
            1: "print(SET)",
            2: "interactive.calculate_distances(SET)",
            3: "interactive.embedding(SET)",
            4: "interactive.plotting(SET)",
            5: "SET = set_collection(tree_set(interactive.add_set())) + SET",  # "interactive.add_set()",
            6: "interactive.get_subset(SET)",
            7: "interactive.exit_program()",
            8: "interactive.usage()",
        }
        return Actions[control]
    except KeyError:
        return "Operation unavailable"
