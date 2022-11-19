__author__ = 'Andrea Rubbi : andrea.rubbi.98@gmail.com'

from tree_set.tree_set import tree_set
from embeddings import PCA_e, tSNE_e
from calculate_distances import hashrf
from interactive_mode import interactive
import tree_emb_parser

 
# ──────────────────────────────────────────────────────────────────────────────
if __name__ =='__main__':
    args = tree_emb_parser.parser()
    SET = tree_set(args.input, args.output, args.distance_matrix, args.metadata)
    print('Your input:')
    print(SET)
    
    if args.hashrf: 
        SET.calculate_distances('hashrf')
    
    if args.pca != None:
        SET.embed('pca', args.pca)
    
    if args.tsne != None:
        SET.embed('tsne', args.tsne)
    
    if args.interactive_mode:
        interactive.usage()
    
    while args.interactive_mode:
        control = int(input('Command: '))
        exec(interactive.interact(control))
    
        
    
    