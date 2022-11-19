import os, sys
# getting the name of the directory
current = os.path.dirname(os.path.realpath(__file__))
 
# Getting the parent directory name
parent = os.path.dirname(current)
 
# adding the parent directory to
# the sys.path.
sys.path.append(parent)

def usage():
    print('TreeEmbedding | Interactive Mode')
    print('⣿⣿⣿⣿⣿⣿⣿⣿⠿⣟⠉⡿⠿⣿⣿⣿⣿⣿⣿⣿ -- Controls --',
          '⣿⣿⣿⣿⡿⣿⢉⢳⠴⣞⠉⡷⢥⡏⡙⡿⢿⣿⣿⣿ 1 --> see tree set status ',
          '⣿⣿⡋⢻⡤⣼⠉⢯⡆⣞⠙⣧⣢⠏⠪⣣⢦⡛⠹⣿ 2 --> calculate distances - hashrf',
          '⣿⣿⠓⢻⣄⣼⠋⢷⡠⡽⠚⣉⣤⡞⢚⢦⢢⠟⠹⣿ 3 --> embed distances - PCA',
          '⣿⣿⢓⢻⡄⡼⠗⢃⣂⡒⠻⣧⣂⡿⠚⣨⢨⡓⠻⣿ 4 --> embed distances - tSNE',
          '⣿⣿⠗⢎⢄⠂⠾⣯⣂⡽⢓⢆⢔⠐⠿⣇⢅⡗⠻⣿ 5 --> exit',
          '⣿⣿⠖⢯⡡⣹⠗⣤⣉⠛⠶⣏⢌⡿⠲⡌⢌⢞⠼⣿ 6 --> see list of controls',
          '⣿⣿⣮⣾⡉⣹⠦⣞⡉⡽⠶⡌⠌⡮⠲⣏⢩⣷⣼⣿',
          '⣿⣿⣿⣿⣿⣿⣤⣞⢉⣳⠥⣏⠍⣧⣵⣿⣿⣿⣿⣿',
          '⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣤⣿⣿⣿⣿⣿⣿⣿⣿⣿', sep = '\n')
    
    return 0

def exit_program():
    print('- Leaving TreeEmbedding -')
    exit()

def ask_dimensions(method, SET):
    dimensions = int(input("Number of components: "))
    if method == 'pca':
        SET.embed('pca', dimensions)
        return 0
    if method == 'tsne':
        SET.embed('tsne', dimensions)
        return 0


def interact(control):
    try:
        Actions = {
            1 : "print(SET)",
            2 : "SET.calculate_distances('hashrf')",
            3 : "interactive.ask_dimensions('pca', SET)",
            4 : "interactive.ask_dimensions('tsne', SET)",
            5 : "interactive.exit_program()",
            6 : "interactive.usage()",
            7 : "",
            8 : "",
            
        }
        return Actions[control]
    except: return "print('Select a key from the list!')"
    