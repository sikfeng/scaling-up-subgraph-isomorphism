package sgi;

import sgi.CountFingerprinter;

import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.exception.*;
import org.openscience.cdk.*;
import org.openscience.cdk.smiles.*;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.fingerprint.*;

import java.io.*;
import java.lang.*;
import java.util.*;

class GraphSig {
  static CountFingerprinter countMolFingerprinter = new CountFingerprinter();
  static IFingerprinter molFingerprinter = new Fingerprinter();
  int[] counted = {6, 18, 14, 15, 8};
  static ArrayList<IAtomContainer> countedGroups = new ArrayList<IAtomContainer>();
  static ArrayList<Pattern> countedGroupPatterns = new ArrayList<Pattern>();
  IBitFingerprint molFingerprint;
  int[] countMolFingerprint;
  int[] count = new int[counted.length];
  int[] countGroup = new int[5];
  static int selSig = 1011; // 111 binary
  // 001 count atom
  // 010 count group
  // 100 fingerprint
  public static void init(int selected) {
    selSig = selected;
    if (((selSig >> 1) & 1) == 1) {
			try {
					SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
					countedGroups.add(sp.parseSmiles("OP=O"));
					//countedGroups.add(sp.parseSmiles("CC"));
			} catch (InvalidSmilesException e) {
					System.err.println(e.getMessage());
					System.exit(1);
			}
			for (int i = 0; i < countedGroups.size(); ++i) {
				countedGroupPatterns.add(Pattern.findSubstructure(countedGroups.get(i)));
			}
		}
  }
  public void update(IAtomContainer mol) {
    if ((selSig & 1) == 1) {
      for (IAtom atom : mol.atoms()) {
        int atomicNum = atom.getAtomicNumber();
        for (int i = 0; i < counted.length; ++i) {
          if (counted[i] == atomicNum) {
            count[i]++;
            break;
          }
        }
      }
    }
    if (((selSig >> 1) & 1) == 1) {
			for (int i = 0; i < countedGroups.size(); ++i) {
				for (int[] mapping : countedGroupPatterns.get(i).matchAll(mol)) {
					countGroup[i]++;
				}
			}
    }
    if (((selSig >> 2) & 1) == 1) {
			try {
				molFingerprint = molFingerprinter.getBitFingerprint(mol);
			} catch (CDKException e) {
					System.err.println(e.getMessage());
					System.exit(1);
			}
    }
    if (((selSig >> 3) & 1) == 1) {
			try {
				countMolFingerprint = countMolFingerprinter.getArrayFingerprint(mol);
			} catch (CDKException e) {
					System.err.println(e.getMessage());
					System.exit(1);
			}
    }
  }
  public boolean compare(GraphSig subgraphSig) {
    if ((selSig & 1) == 1) {
			for (int i = 0; i < counted.length; ++i) {
				if (subgraphSig.count[i] > this.count[i]) {
					return false;
				}
			}
		}
    if (((selSig >> 1) & 1) == 1) {
			for (int i = 0; i < countedGroups.size(); ++i) {
				if (subgraphSig.countGroup[i] > this.countGroup[i]) {
					return false;
				}
			}
    }
    if (((selSig >> 2) & 1) == 1) {
			for (int i = 0; i < Math.min(this.molFingerprint.size(), subgraphSig.molFingerprint.size()); ++i) {
			 if (!this.molFingerprint.get(i) && subgraphSig.molFingerprint.get(i)) {
				 return false;
			 }
			}
    }
    if (((selSig >> 3) & 1) == 1) {
      for (int i = 0; i < Math.min(this.countMolFingerprint.length, subgraphSig.countMolFingerprint.length); ++i) {
        if (this.countMolFingerprint[i] < subgraphSig.countMolFingerprint[i]) {
          return false;
        }
      }
    }
		return true;
  }
}

