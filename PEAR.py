__author__ = "Andrea Rubbi : andrea.rubbi.98@gmail.com"

import pear_ebi.tree_emb_parser
from pear_ebi.calculate_distances import hashrf
from pear_ebi.embeddings import PCA_e, tSNE_e
from pear_ebi.interactive_mode import interactive
from pear_ebi.tree_set import set_collection, tree_set

# ──────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    args = pear_ebi.tree_emb_parser.parser()
    if type(args.input) == type(list()):
        SET = set_collection(args.input, args.output, args.distance_matrix, args.metadata)
    else:
        SET = tree_set(args.input, args.output, args.distance_matrix, args.metadata)

    if args.verbose or args.interactive_mode:
        print("Your input:")
        print(SET)

    if args.hashrf:
        SET.calculate_distances("hashrf")

    if args.hashrf_weighted:
        SET.calculate_distances("hashrf_weighted")

    if args.days_RF:
        SET.calculate_distances("days_RF")

    if args.pca != None:
        SET.embed("pca", args.pca)

    if args.tsne != None:
        SET.embed("tsne", args.tsne)

    if args.plot != None:
        method = "tsne" if args.tsne is not None else "pca"
        if args.plot == 2:
            plot = SET.plot_2D(method, save=True)
        if args.plot == 3:
            plot = SET.plot_3D(method, save=True)

    if args.subset != None:
        SET.get_subset(args.subset)

    if args.interactive_mode:
        interactive.usage()

    while args.interactive_mode:
        control = input("Command: ")
        try:
            control = int(control)
        except:
            pass
        exec(interactive.interact(control))
