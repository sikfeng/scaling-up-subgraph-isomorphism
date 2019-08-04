package sgi;

import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.exception.*;
import org.openscience.cdk.*;
import org.openscience.cdk.smiles.*;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.fingerprint.*;
import org.openscience.cdk.similarity.Tanimoto;

import java.io.*;
import java.lang.*;
import java.util.*;

class GraphSig {
  static IFingerprinter molFingerprinter = new Fingerprinter();
  IBitFingerprint molFingerprint;
  static int selSig = 0;
  public static void init(int selected) {
    selSig = selected;
  }
  public void update(IAtomContainer mol) {
    try {
      molFingerprint = molFingerprinter.getBitFingerprint(mol);
    } catch (CDKException e) {
        System.err.println(e.getMessage());
        System.exit(1);
    }
  }
  public boolean compare(GraphSig subgraphSig) {
    double similarity = Tanimoto.calculate(this.molFingerprint, subgraphSig.molFingerprint);
    if (similarity < (double) selSig/100) {
      return false;
    }
		return true;
  }
}

