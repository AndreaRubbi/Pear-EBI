from .. import tree_set

class set_collection:
    def __init__(self, set_collection = list(), file = 'Set_collection'):
        for element in set_collection:
            assert type(element)=='tree_set', 'Every element in a set_collection must be a tree_set'
        
        self.set_collection = set_collection
        self.file = file
        
        if len(self.set_collection) > 0:
            with open(self.file, 'w') as trees:
                for set in self.set_collection:
                    trees.write(set.file.read())
                trees.close()
    
    def __add__(self, other):
        if type(other)=='list':
            for element in other:
                assert type(element)=='tree_set', 'Cannot concatenate set_collection with {wrong}'.format(wrong=type(element))
            self.set_collection.expand(other)
            
            with open(self.file, 'a+') as trees:
                for set in other:
                    trees.write(set.file.read())
                trees.close()

        else: 
            assert type(other)=='tree_set', 'Cannot concatenate set_collection with {wrong}'.format(wrong=type(other))
            self.set_collection.append(other)
            with open(self.file, 'a+') as trees:
                trees.write(other.file.read())
                trees.close()
                
    def __str__(self):
        pass
        
    
    def concatenate(self, other):
        if type(other)=='list':
            for element in other:
                assert type(element)=='tree_set', 'Cannot concatenate set_collection with {wrong}'.format(wrong=type(element))
            self.set_collection.expand(other)
            
            with open(self.file, 'a+') as trees:
                for set in other:
                    trees.write(set.file.read())
                trees.close()

        else: 
            assert type(other)=='tree_set', 'Cannot concatenate set_collection with {wrong}'.format(wrong=type(other))
            self.set_collection.append(other)
            with open(self.file, 'a+') as trees:
                trees.write(other.file.read())
                trees.close()
            
            
            
        
        

         