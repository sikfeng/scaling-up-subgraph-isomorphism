import igraph 
import pickle 
import time 
import sys

def unpickle(fpath): 
    return pickle.load(open(fpath, 'rb')) 

startExtractFile = time.time() 
graphs = unpickle("graphs" + str(sys.argv[1]) + ".pkl") # list of graphs 
subgraphs = unpickle("subgraphs" + str(sys.argv[1]) + ".pkl") # list of subgraphs 
print("Extracted structures in {:.3f} seconds".format(time.time() - startExtractFile)) 
print("{} graphs and {} subgraphs".format(len(graphs), len(subgraphs)))

# direct SGI
print("\ndirect")
startSGITime = time.time() 
numSGI = 0
def atomToInt(atoms): 
    # because subisomorphic only accepts integers for color 
    intColor = [] 
    for atom in atoms: 
        num = ord(atom[0]) - 65 # ord('A') == 65 
        if len(atom) == 2: 
            num *= 100 
            num += ord(atom[1]) - 97 # ord('a') == 97 
        intColor.append(num) 
    return intColor 
isomorphicDirect = [] 
for i in range(len(graphs)): 
    g = graphs[i] 
    for j in range(len(subgraphs)): 
        sg = subgraphs[j] 
        numSGI += 1
        sgiso = g.subisomorphic_vf2(sg, color1=atomToInt(g.vs["atom"]), color2=atomToInt(sg.vs["atom"]), edge_color1=g.es['weight'], edge_color2=sg.es['weight']) 
        if sgiso: 
            isomorphicDirect.append((i, j,)) 
print("Finished direct SGI in {:.3f} seconds".format(time.time() - startSGITime))
print("Found {} substructures, SGI run {} times ".format(len(isomorphicDirect), numSGI))


# SGI with graph signature
print("\ngraph signature")
startSGITime = time.time()
class graphSig: 
    def __init__(self, g): 
        self.molCount = {'C': 0, 'S': 0, 'P': 0, 'O': 0, 'Si': 0} 
        self.color = [] 
        atoms = g.vs['atom'] 
        for atom in atoms: 
            if atom in self.molCount: 
                self.molCount[atom] += 1 
            # convert atom to integer 
            # because subisomorphic only accepts integers for color 
            num = ord(atom[0]) - 65 # ord('A') == 65 
            if len(atom) == 2: 
                num *= 100 
                num += ord(atom[1]) - 97 # ord('a') == 97 
            self.color.append(num) 
    def compare(self, subgraphSig): 
        for atom in self.molCount: 
            if self.molCount[atom] < subgraphSig.molCount[atom]: 
                return False 
        return True 
isomorphicSig = [] 
gSig = [] 
subgSig = [] 
numSGI = 0
for i in range(len(graphs)): 
    g = graphs[i] 
    gSig.append(graphSig(g)) 
    for j in range(len(subgraphs)): 
        sg = subgraphs[j] 
        if i == 0: 
            subgSig.append(graphSig(sg)) 
        sigCheck = gSig[i].compare(subgSig[j]) 
        if sigCheck:
            numSGI += 1
            sgiso = g.subisomorphic_vf2(sg, color1=gSig[i].color, color2=subgSig[j].color, edge_color1=g.es['weight'], edge_color2=sg.es['weight']) 
            if sgiso: 
                isomorphicSig.append((i, j,)) 
print("Finished SGI with graph signatures in {:.3f} seconds".format(time.time() - startSGITime)) 
print("Found {} substructures, SGI run {} times ".format(len(isomorphicSig), numSGI))
