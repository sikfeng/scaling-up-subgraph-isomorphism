package sgi;

import sgi.CountFingerprinter;

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

public class GraphSig {

  static CountFingerprinter molFingerprinter = new CountFingerprinter();
  int[] molFingerprint;
  static final int[] countedAtoms = {6, 8, 14, 15, 18};
  int[] countAtoms = new int[countedAtoms.length];

  public static void init() {
    // currently nothing needed
  }

  public void update(IAtomContainer mol) {
    try {
      molFingerprint = molFingerprinter.getArrayFingerprint(mol);
    } catch (CDKException e) {
        e.printStackTrace();
        System.exit(1);
    }
    Iterable<IAtom> atoms = mol.atoms();

    for (IAtom atom : atoms) {
      int atomicNum = atom.getAtomicNumber();
      for (int i = 0; i < countedAtoms.length; ++i) {
        if (countedAtoms[i] == atomicNum) {
          countAtoms[i]++;
          break;
        }
      }
    }
  }

  public int compare(GraphSig subgraphSig) {
    for (int i = 0; i < countedAtoms.length; ++i) {
      if (subgraphSig.countAtoms[i] > this.countAtoms[i]) {
        return -1;
      }
    }
    
    for (int i = 0; i < Math.min(this.molFingerprint.length, subgraphSig.molFingerprint.length); ++i) {
      if (this.molFingerprint[i] < subgraphSig.molFingerprint[i]) {
        return -1;
      }
    }
    
    return 0;
  }
}

