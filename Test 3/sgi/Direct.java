package sgi;

import java.io.*;
import java.lang.*;
import java.util.*;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.isomorphism.Pattern;

public class Direct {
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
      // direct SGI
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
    } catch (Exception e) { // unreported exception errors and stuff
      System.out.println("Error: " + e.getMessage());
      System.exit(1);
    } 
  }
}
