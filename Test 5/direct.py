import sys
import time
from rdkit import Chem
from multiprocessing import Pool


def main():
    start_extract_file = time.time()
    graphsFile = Chem.SDMolSupplier('sdf/candidates/' + str(sys.argv[1]) + '/0.sdf',sanitize=False)
    subgraphsFile = Chem.SDMolSupplier('sdf/substructures_all/0.sdf',sanitize=False)
    gs = [x for x in graphsFile]
    global sgs
    sgs = [x for x in subgraphsFile]
    print("Extracted structures in {:.2f} seconds".format(time.time() - start_extract_file))

    start_sgi = time.time()
    numSGIRun = 0
    numSubIso = 0
    with Pool(processes=8) as pool:
        results = list(pool.imap_unordered(find_subgraphs, gs, chunksize=250))
    for subIso, SGI in results:
        numSGIRun += SGI
        numSubIso += subIso
    print("Finished SGI in {:.2f} seconds".format(time.time() - start_sgi))
    print("Found {} substructures, SGI run {} times".format(numSubIso, numSGIRun))

def find_subgraphs(graph):
    numSubIso = 0
    numSGIRun = 0
    for i in range(len(sgs)):
        numSGIRun += 1
        subgraph = sgs[i]
        if graph.HasSubstructMatch(subgraph):
            numSubIso += 1
    return (numSubIso, numSGIRun)

main()
