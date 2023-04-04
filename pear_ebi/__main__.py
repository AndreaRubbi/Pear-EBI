import os
import sys


def main():
    __author__ = "Andrea Rubbi : andrea.rubbi.98@gmail.com"

    import os
    import re
    import sys
    from collections import defaultdict
    from glob import glob

    import numpy as np
    import pandas as pd
    import toml
    from rich import print

    # getting the name of the directory
    current = os.path.dirname(os.path.realpath(__file__))

    # Getting the parent directory name
    parent = os.path.dirname(current)

    # adding the parent directory to
    # the sys.path.
    sys.path.append(parent)

    import pear_ebi.tree_emb_parser
    from pear_ebi.calculate_distances import hashrf
    from pear_ebi.embeddings import PCA_e, tSNE_e
    from pear_ebi.interactive_mode import interactive
    from pear_ebi.tree_set import set_collection, tree_set

    try:
        # Retrieves args from parser
        args = pear_ebi.tree_emb_parser.parser()

        if args.version:
            print(f"[blue]PEAR v{pear_ebi.__version__}")

        # ─── Interactive Mode ─────────────────────────────────────────────────
        if args.interactive_mode:
            file = args.input
            if not file:
                while True:
                    try:
                        print("[bright_magenta]Specify file with tree set\n")
                        file.append(input("File: "))
                        break
                    except FileNotFoundError:
                        print("[red] File not found - try again")
                    except ValueError:
                        print("[red] Invalid file specification")

            distance_matrix = args.distance_matrix
            output_file = args.output
            metadata = args.metadata

            if len(file) == 1:
                SET = tree_set(
                    file[0],
                    output_file=output_file,
                    distance_matrix=distance_matrix,
                    metadata=metadata,
                )

            else:
                SET = set_collection(
                    collection=file,
                    output_file=output_file,
                    distance_matrix=distance_matrix,
                    metadata=metadata,
                )

            # shows set specifics
            print("[bright_magenta]Your input:")
            print(SET)

            interactive.usage()

        while args.interactive_mode:
            control = ""
            try:
                control = input("Command: ")
                control = int(control)
            except ValueError:
                pass
            except KeyboardInterrupt:
                print("[orange1]\n- Leaving PEAR -")
                exit()
            try:
                exec(interactive.interact(control))
            except KeyboardInterrupt:
                print("[red] Interrupted")
                pass

        # ─── Normal Pear ──────────────────────────────────────────────────────
        else:
            # if config file is specified,
            # pear looks for it and tries to
            # load the configurations
            config = defaultdict(lambda: None)
            if args.config is not None:
                try:
                    config_file = args.config
                    config = toml.load(config_file)
                except FileNotFoundError:
                    print("[red]File not found")
                    exit()
                except TypeError:
                    print("[red]Config arguments must be a filename or filepath")
                    exit()
            else:
                # Automatically check whether there is
                # a pear.toml file in directory
                try:
                    config = toml.load("./pear.toml")
                # If pear.toml not found, then continue
                except FileNotFoundError:
                    pass

            if "config" in globals() or "config" in locals():  # vars(__builtins__):
                config = defaultdict(lambda: None, config)
                # print(config)

            # files is a list that is populated from
            # the config file and the parser,
            # it starts by iterating all the files
            # in the config file at key files,
            # then iterating all the files in all
            # the directories at key dir,
            # then iterating all the files
            # in arg.inputs

            # ─── Define Tree Set Or Collection ────────────────────────────
            files = list()
            if config["dir"] is not None:
                config["dir"] = defaultdict(lambda: None, config["dir"])
                pattern = (
                    config["dir"]["pattern"]
                    if config["dir"]["pattern"] is not None
                    else "*"
                )
                pattern = args.pattern if args.pattern is not None else pattern
            else:
                pattern = args.pattern if args.pattern is not None else "*"

            if config["trees"] is not None:
                for t in config["trees"].values():
                    files.append(t)

            if config["dir"] is not None:
                for k, d in config["dir"].items():
                    if k[:3] != "dir":
                        continue
                    else:
                        print(
                            f"[white]Looking into directory [purple]{d} [white]- pattern: [purple]{pattern}"
                        )
                    try:
                        glob_pattern = os.path.join(d, pattern)
                        files.extend(glob(glob_pattern))
                    except FileNotFoundError:
                        print("[red]File or directory not found")
                    except ValueError:
                        print("[red]Invalid type")
            if args.dir is not None:
                try:
                    pattern = args.pattern if args.pattern is not None else "*"
                    glob_pattern = os.path.join(args.dir, pattern)
                    files.extend(glob(glob_pattern))
                except FileNotFoundError:
                    print("[red]File or directory not found")
                except ValueError:
                    print("[red]Invalid type")

            files.extend(args.input)

            # pear will simply close if no file
            # is given!
            if not files:
                print("[red]No files specified[white] (see --help for instructions)")
                print("[orange1]- Leaving PEAR -")
                exit()

            # if the key collection is present, then
            # pear parses its keys:
            # output_file, distance_matrix,
            # and metadata
            if config["collection"] is not None:
                config["collection"] = defaultdict(lambda: None, config["collection"])
            else:
                config["collection"] = defaultdict(lambda: None)

            distance_matrix = config["collection"]["distance_matrix"]
            output_file = config["collection"]["output_file"]
            metadata = config["collection"]["metadata"]

            # if the same args are present in
            # the pear call, they are
            # substituted by these
            distance_matrix = (
                args.distance_matrix
                if args.distance_matrix is not None
                else distance_matrix
            )
            output_file = args.output if args.output is not None else output_file
            metadata = args.metadata if args.metadata is not None else metadata

            # we can now define our set_collection
            if len(files) == 1:
                SET = tree_set(
                    files[0],
                    output_file=output_file,
                    distance_matrix=distance_matrix,
                    metadata=metadata,
                )

            else:
                SET = set_collection(
                    collection=list(map(lambda f: tree_set(f), files)),
                    output_file=output_file,
                    distance_matrix=distance_matrix,
                    metadata=metadata,
                )

            # shows set specifics
            print("[bright_magenta]Your input:")
            print(SET)

            # same as above, we try to read method
            # value from config and we overwrite it
            # with the value given in terminal
            # ─── Compute Distances ────────────────────────────────────────
            method = config["method"] if args.method is None else args.method
            if method is not None:
                SET.calculate_distances(method)

            # embedding is optional, so we check
            # if any value is given in the config file
            # every parameter is a list
            # because we allow for multiple embeddings
            # ─── Compute Embeddings ───────────────────────────────────────
            if config["embedding"] is not None:
                config["embedding"] = defaultdict(lambda: None, config["embedding"])
                method_embedding = (
                    config["embedding"]["method"]
                    if config["embedding"]["method"] is not None
                    else "pca"
                )
                dimensions = (
                    config["embedding"]["dimensions"]
                    if config["embedding"]["dimensions"] is not None
                    else 2
                )
                quality = (
                    config["embedding"]["quality"]
                    if config["embedding"]["quality"] is not None
                    else False
                )
                report = (
                    config["embedding"]["report"]
                    if config["embedding"]["report"] is not None
                    else False
                )

            # tries to overwrite
            # these parameters
            if config["embedding"] is None:
                quality = args.quality
            if args.quality:
                quality = True

            if config["embedding"] is None:
                report = args.report
                if report:
                    quality = True
            if args.report:
                quality = True
                report = True

            if args.pca is not None:
                method_embedding = "pca"
                dimensions = args.pca
            elif config["embedding"] is None and args.tsne is None:
                method_embedding = None
                dimensions = None

            if args.tsne is not None:
                method_embedding = "tsne"
                dimensions = args.tsne
            elif config["embedding"] is None and args.pca is None:
                method_embedding = None
                dimensions = None

            # parse values and returns embeddings
            if method_embedding is not None:
                SET.embed(
                    method=method_embedding,
                    dimensions=dimensions,
                    quality=quality,
                    report=report,
                )
                if report:
                    SET.emb_quality.report()

            # ─── Highlights ───────────────────────────────────────
            if config["highlight"] is not None and config["trees"] is not None:
                config["highlight"] = defaultdict(lambda: None, config["highlight"])
                highlight = list()
                tree_files_config = list(
                    map(
                        lambda f: os.path.splitext(os.path.basename(f))[0],
                        config["trees"].values(),
                    )
                )

                n_key_tree_config = list(
                    map(lambda k: re.findall(r"\d+", k)[0], config["trees"].keys())
                )

                for set_trees in SET.data.keys():
                    if config["highlight"][set_trees] is not None:
                        values = np.array(
                            [0 for i in range(SET.data[set_trees]["n_trees"])]
                        )
                        values[config["highlight"][set_trees]] = 1
                        highlight.extend(values.tolist())

                    elif set_trees in tree_files_config:
                        idx = tree_files_config.index(set_trees)
                        if (
                            config["highlight"][f"trace{n_key_tree_config[idx]}"]
                            is not None
                        ):
                            values = np.array(
                                [0 for i in range(SET.data[set_trees]["n_trees"])]
                            )
                            values[
                                config["highlight"][f"trace{n_key_tree_config[idx]}"]
                            ] = 1
                            highlight.extend(values.tolist())

                        else:
                            highlight.extend(
                                [0 for i in range(SET.data[set_trees]["n_trees"])]
                            )

                SET.metadata["highlight"] = highlight

            # ─── Plot Embeddings ──────────────────────────────────────────
            if method_embedding is not None:
                if config["plot"] is not None:
                    config["plot"] = (
                        defaultdict(lambda: None, config["plot"])
                        if config["plot"] is not None
                        else defaultdict(lambda: None)
                    )
                    name_plot = config["plot"]["name"]
                    plot_meta = (
                        config["plot"]["plot_meta"]
                        if config["plot"]["meta"] is not None
                        else "SET-ID"
                    )
                    plot_set = config["plot"]["plot_set"]
                    select = (
                        config["plot"]["select"]
                        if config["plot"]["select"] is not None
                        else False
                    )
                    same_scale = (
                        config["plot"]["same_scale"]
                        if config["plot"]["same_scale"] is not None
                        else False
                    )

                    show = (
                        config["plot"]["show"]
                        if config["plot"]["show"] is not None
                        else False
                    )

                    if dimensions > 2:
                        name_plot = (
                            name_plot + "3D"
                            if name_plot is not None
                            else f"{method_embedding.upper()}_3D"
                        )
                        fig = SET.plot_3D(
                            method_embedding,
                            name_plot=name_plot,
                            plot_meta=plot_meta,
                            plot_set=plot_set,
                            select=select,
                            same_scale=same_scale,
                            save=True,
                        )

                        if show:
                            fig.show()

                    name_plot = (
                        name_plot + "2D"
                        if name_plot is not None
                        else f"{method_embedding.upper()}_2D"
                    )
                    fig = SET.plot_2D(
                        method_embedding,
                        name_plot=name_plot,
                        plot_meta=plot_meta,
                        plot_set=plot_set,
                        select=select,
                        same_scale=same_scale,
                        save=True,
                    )

                    if show:
                        fig.show()

                else:
                    if dimensions > 2:
                        name_plot = f"{method_embedding.upper()}_3D"
                        fig = SET.plot_3D(
                            method_embedding,
                            name_plot=name_plot,
                            save=True,
                        )

                        if args.plot:
                            fig.show()

                        fig.show()
                    name_plot = f"{method_embedding.upper()}_2D"
                    fig = SET.plot_2D(
                        method_embedding,
                        name_plot=name_plot,
                        save=True,
                    )

                    if args.plot:
                        fig.show()

            # ─── Get Subset ───────────────────────────────────────────────
            if args.subset != None:
                fig2, fig3 = SET.get_subset(args.subset)
                fig2.show()
                fig3.show()
    except KeyboardInterrupt:
        print("[orange1]\n- Leaving PEAR -")

    print("[orange1]\n- Leaving PEAR -")


# ──────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    sys.exit(main())
