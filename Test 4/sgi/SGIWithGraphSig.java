package sgi;

import sgi.GraphSig;
import java.io.*;
import java.lang.*;
import java.util.*;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.isomorphism.Pattern;

public class SGIWithGraphSig {
  public static void main(String[] args) {
    try {
      double startExtractTime = System.currentTimeMillis();
      File graphsFile = new File("graphs" + args[0] + "/0.sdf");
      File subgraphsFile = new File("subgraphs" + args[0] + "/0.sdf");
      IteratingSDFReader graphReader = new IteratingSDFReader(new FileInputStream(graphsFile), DefaultChemObjectBuilder.getInstance());
      IteratingSDFReader subgraphReader = new IteratingSDFReader(new FileInputStream(subgraphsFile), DefaultChemObjectBuilder.getInstance());
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
      // SGI with graph signatures
      GraphSig.init(Integer.parseInt(args[1]));
      double startSGITimeSig = System.currentTimeMillis();;
      int numSubIsoSig = 0;
      int numSGIRunSig = 0;

      ArrayList<GraphSig> graphSigs = new ArrayList<GraphSig>();
      ArrayList<GraphSig> subgraphSigs = new ArrayList<GraphSig>();

      for (int i = 0; i < graphs.size(); ++i) {
        GraphSig graphSig = new GraphSig();
        graphSig.update(graphs.get(i));
        graphSigs.add(graphSig);
      }
      for (int j = 0; j < subgraphs.size(); ++j) {
        GraphSig subgraphSig = new GraphSig();
        subgraphSig.update(subgraphs.get(j));
        subgraphSigs.add(subgraphSig);
      }
      for (int i = 0; i < graphs.size(); ++i) {
        GraphSig graphSig = graphSigs.get(i);
        for (int j = 0; j < subgraphs.size(); ++j) {
          GraphSig subgraphSig = subgraphSigs.get(j);
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
