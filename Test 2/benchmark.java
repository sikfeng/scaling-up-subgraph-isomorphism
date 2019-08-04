import java.io.*;
import java.lang.*;
import java.util.*;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.isomorphism.Pattern;

class GraphSig {
  int[] counted = {6, 18, 14, 15, 8};
  int[] count = new int[counted.length];
  public void update(IAtomContainer mol) {
    Iterable<IAtom> atoms = mol.atoms();
    for (IAtom atom : atoms) {
      int atomicNum = atom.getAtomicNumber();
      for (int i = 0; i < counted.length; ++i) {
        if (counted[i] == atomicNum) {
          count[i]++;
        }
      }
    }
  }
  public boolean compare(GraphSig subgraphSig) {
    for (int i = 0; i < counted.length; ++i) {
      if (subgraphSig.count[i] > this.count[i]) {
        return false;
      }
    }
    return true;
  }
}

public class benchmark {
  public static void main(String[] args) {
    try {
      double startExtractTime = System.currentTimeMillis();
      File graphsFile = new File("graphs1/0.sdf");
      File subgraphsFile = new File("subgraphs1/0.sdf");
      IteratingSDFReader graphReader = new IteratingSDFReader(new FileInputStream(subgraphsFile), DefaultChemObjectBuilder.getInstance());
      IteratingSDFReader subgraphReader = new IteratingSDFReader(new FileInputStream(graphsFile), DefaultChemObjectBuilder.getInstance());
      ArrayList<IAtomContainer> graphs = new ArrayList<IAtomContainer>();
      ArrayList<IAtomContainer> subgraphs = new ArrayList<IAtomContainer>();
      while (graphReader.hasNext()) {
        IAtomContainer graph = graphReader.next();
        graphs.add(graph);
      }
      while (subgraphReader.hasNext()) {
        IAtomContainer subgraph = subgraphReader.next();
        subgraphs.add(subgraph);
      }
      System.out.println("Extracted structures in " + ((System.currentTimeMillis() - startExtractTime) / 1000) + " seconds");
      System.out.println(graphs.size() + " graphs, " + subgraphs.size() + " subgraphs");

      // direct SGI
      System.out.println();
      System.out.println("direct");
      double startSGITimeDirect = System.currentTimeMillis();;
      int numSubIsoDirect = 0;
      int numSGIRunDirect = 0;
      for (int i = 0; i < graphs.size(); ++i) {
        for (int j = 0; j < subgraphs.size(); ++j) {
          Pattern subgraphPattern = Pattern.findSubstructure(subgraphs.get(j));
          boolean isSubgraph = subgraphPattern.matches(graphs.get(i));
          ++numSGIRunDirect;
          if (isSubgraph) {
            ++numSubIsoDirect;
          }
        }
      }
      System.out.println("Finished SGI in " + ((System.currentTimeMillis() - startSGITimeDirect) / 1000) + " seconds");
      System.out.println("Found " + numSubIsoDirect + " substructures, SGI run " + numSGIRunDirect + " times");

      // SGI with graph signatures
      System.out.println();
      System.out.println("graph signatures");
      double startSGITimeSig = System.currentTimeMillis();;
      int numSubIsoSig = 0;
      int numSGIRunSig = 0;
      ArrayList<GraphSig> graphSigs = new ArrayList<GraphSig>();
      ArrayList<GraphSig> subgraphSigs = new ArrayList<GraphSig>();
      for (int i = 0; i < graphs.size(); ++i) {
        GraphSig graphSig = new GraphSig();
        graphSig.update(graphs.get(i));
        graphSigs.add(graphSig);
        for (int j = 0; j < subgraphs.size(); ++j) {
          if (i == 0) { // only run on first time when j is 0, then memoize
            GraphSig subgraphSig = new GraphSig();
            subgraphSig.update(subgraphs.get(j));
            subgraphSigs.add(subgraphSig);
          }
          if (graphSig.compare(subgraphSigs.get(j))) {
            ++numSGIRunSig;
            Pattern subgraphPattern = Pattern.findSubstructure(subgraphs.get(j));
            boolean isSubgraph = subgraphPattern.matches(graphs.get(i));
            if (isSubgraph) {
              ++numSubIsoSig;
            }
          }
        }
      }
      System.out.println("Finished SGI in " + ((System.currentTimeMillis() - startSGITimeSig) / 1000) + " seconds");
      System.out.println("Found " + numSubIsoSig + " substructures, SGI run " + numSGIRunSig + " times");
    } catch (Exception e) { // unreported exception errors and stuff
      System.out.println("Error: " + e.getMessage());
      System.exit(1);
    } 
  }
}
