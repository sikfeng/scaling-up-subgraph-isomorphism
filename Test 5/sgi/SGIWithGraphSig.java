package sgi;

import sgi.CountFingerprinter;
import sgi.GraphSig;

import java.io.*;
import java.lang.*;
import java.lang.Math.*;
import java.util.*;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.exception.*;
import org.openscience.cdk.*;
import org.openscience.cdk.smiles.*;
import org.openscience.cdk.fingerprint.*;
import org.openscience.cdk.tools.*;
import org.openscience.cdk.tools.manipulator.*;
import org.openscience.cdk.config.*;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.isomorphism.Pattern;

public class SGIWithGraphSig {
  public static void main(String[] args) {
    try {
      double startExtractTime = System.currentTimeMillis();
      File graphsFile = new File("sdf/candidates/" + args[0] + "/0.sdf");
      File subgraphsFile = new File("sdf/substructures_all/0.sdf");
      IteratingSDFReader graphReader = new IteratingSDFReader(new FileInputStream(graphsFile), SilentChemObjectBuilder.getInstance());
      IteratingSDFReader subgraphReader = new IteratingSDFReader(new FileInputStream(subgraphsFile), SilentChemObjectBuilder.getInstance());
      ArrayList<IAtomContainer> subgraphs = new ArrayList<IAtomContainer>();
      ArrayList<IAtomContainer> graphs = new ArrayList<IAtomContainer>();
      while (subgraphReader.hasNext()) {
        subgraphs.add(subgraphReader.next());
      }
      while (graphReader.hasNext()) {
        graphs.add(graphReader.next());
      }
      System.out.println("Extracted structures in " + ((System.currentTimeMillis() - startExtractTime) / 1000) + " seconds");
      
			double startSigTime = System.currentTimeMillis();
      GraphSig.init();
      ArrayList<GraphSig> subgraphSigs = new ArrayList<GraphSig>();
      ArrayList<GraphSig> graphSigs = new ArrayList<GraphSig>();

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
      System.out.println("Generated signatures in " + ((System.currentTimeMillis() - startExtractTime) / 1000) + " seconds");
      
      double startSGITime = System.currentTimeMillis();
      int numSGIRun = 0;
     	int numSubIso = 0; 
			for (int i = 0; i < graphs.size(); ++i) {
        IAtomContainer graph = graphs.get(i);
        GraphSig graphSig = graphSigs.get(i);
        for (int j = 0; j < subgraphs.size(); ++j) {
          if (graphSig.compare(subgraphSigs.get(j)) == 0) {
            ++numSGIRun;
            Pattern subgraphPattern = Pattern.findSubstructure(subgraphs.get(j));
            boolean isSubgraph = subgraphPattern.matches(graph);
            if (isSubgraph) {
              ++numSubIso;
            }
          }
        }
      }
      System.out.println("Finished SGI in " + ((System.currentTimeMillis() - startSGITime) / 1000) + " seconds");
      System.out.println("Found " + numSubIso + " substructures, SGI run " + numSGIRun
                      + " times");
      System.out.println(graphs.size() + " num of graphs");
      System.out.println(subgraphs.size() + " num of subgraphs");
    } catch (Exception e) { // unreported exception errors and stuff
      System.out.println("Error: " + e);
      System.exit(1);
    } 
  }
}
