package sgi;

import java.io.*;
import java.lang.*;
import java.util.*;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.exception.*;
import org.openscience.cdk.*;
import org.openscience.cdk.smiles.*;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.fingerprint.*;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.AbstractMap.SimpleImmutableEntry;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

public class CountFingerprinter extends AbstractFingerprinter implements IFingerprinter {

    /** Throw an exception if too many paths (per atom) are generated. */
    private final static int                 DEFAULT_PATH_LIMIT   = 42000;

    /** The default length of created fingerprints. */
    public final static int                  DEFAULT_SIZE         = 1024;
    /** The default search depth used to create the fingerprints. */
    public final static int                  DEFAULT_SEARCH_DEPTH = 7;

    private int                              size;
    private int                              searchDepth;
    private int                              pathLimit = DEFAULT_PATH_LIMIT;

    private boolean                          hashPseudoAtoms = false;

    static int                               debugCounter         = 0;


    private static ILoggingTool              logger               = LoggingToolFactory
                                                                          .createLoggingTool(Fingerprinter.class);



    /**
     * Creates a fingerprint generator of length <code>DEFAULT_SIZE</code>
     * and with a search depth of <code>DEFAULT_SEARCH_DEPTH</code>.
     */
    public CountFingerprinter() {
        this(DEFAULT_SIZE, DEFAULT_SEARCH_DEPTH);
    }

    public CountFingerprinter(int size) {
        this(size, DEFAULT_SEARCH_DEPTH);
    }

    /**
     * Constructs a fingerprint generator that creates fingerprints of
     * the given size, using a generation algorithm with the given search
     * depth.
     *
     * @param  size        The desired size of the fingerprint
     * @param  searchDepth The desired depth of search (number of bonds)
     */
    public CountFingerprinter(int size, int searchDepth) {
        this.size = size;
        this.searchDepth = searchDepth;

    }

    @Override
    protected List<Map.Entry<String, String>> getParameters() {
        return Arrays.<Map.Entry<String,String>>asList(
            new SimpleImmutableEntry<>("searchDepth", Integer.toString(searchDepth)),
            new SimpleImmutableEntry<>("pathLimit", Integer.toString(pathLimit)),
            new SimpleImmutableEntry<>("hashPseudoAtoms", Boolean.toString(hashPseudoAtoms))
        );
    }

    /**
     * Generates a fingerprint of the default size for the given AtomContainer.
     *
     * @param container The AtomContainer for which a Fingerprint is generated
     * @param ringFinder An instance of
     *                   {@link org.openscience.cdk.ringsearch.AllRingsFinder}
     * @exception CDKException if there is a timeout in ring or aromaticity
     *                         perception
     * @return A {@link BitSet} representing the fingerprint
     */
    public IBitFingerprint getBitFingerprint(IAtomContainer container, AllRingsFinder ringFinder) throws CDKException {
        logger.debug("Entering Fingerprinter");
        logger.debug("Starting Aromaticity Detection");
        long before = System.currentTimeMillis();
        if (!hasPseudoAtom(container.atoms())) {
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(container);
            Aromaticity.cdkLegacy().apply(container);
        }
        long after = System.currentTimeMillis();
        logger.debug("time for aromaticity calculation: " + (after - before) + " milliseconds");
        logger.debug("Finished Aromaticity Detection");
        BitSet bitSet = new BitSet(size);
        //encodePaths(container, searchDepth, bitSet, size);
        return new BitSetFingerprint(bitSet);
    }
    //public IBitFingerprint getBitFingerprint(IAtomContainer container, AllRingsFinder ringFinder) throws CDKException {
    public int[] getArrayFingerprint(IAtomContainer container, AllRingsFinder ringFinder) throws CDKException {
        logger.debug("Entering Fingerprinter");
        logger.debug("Starting Aromaticity Detection");
        long before = System.currentTimeMillis();
        if (!hasPseudoAtom(container.atoms())) {
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(container);
            Aromaticity.cdkLegacy().apply(container);
        }
        long after = System.currentTimeMillis();
        logger.debug("time for aromaticity calculation: " + (after - before) + " milliseconds");
        logger.debug("Finished Aromaticity Detection");
        //BitSet bitSet = new BitSet(size);
        int[] arrayList = new int[size];
        //encodePaths(container, searchDepth, bitSet, size);
        encodePaths(container, searchDepth, arrayList, size);
        //return new BitSetFingerprint(bitSet);
        return arrayList;
    }

    /**
     * Generates a fingerprint of the default size for the given AtomContainer.
     *
     *@param container The AtomContainer for which a Fingerprint is generated
     */
    @Override
    public IBitFingerprint getBitFingerprint(IAtomContainer container) throws CDKException {
        return getBitFingerprint(container, null);
    }
    public int[] getArrayFingerprint(IAtomContainer container) throws CDKException {
        return getArrayFingerprint(container, null);
    }

    /** {@inheritDoc} */
    @Override
    public Map<String, Integer> getRawFingerprint(IAtomContainer iAtomContainer) throws CDKException {
        throw new UnsupportedOperationException();
    }

    private IBond findBond(List<IBond> bonds, IAtom beg, IAtom end) {
        for (IBond bond : bonds)
            if (bond.contains(beg) && bond.contains(end))
                return bond;
        return null;
    }

    private String encodePath(IAtomContainer mol, Map<IAtom, List<IBond>> cache, List<IAtom> path, StringBuilder buffer) {
        buffer.setLength(0);
        IAtom prev = path.get(0);
        buffer.append(getAtomSymbol(prev));
        for (int i = 1; i < path.size(); i++) {
            final IAtom next  = path.get(i);
            List<IBond> bonds = cache.get(prev);

            if (bonds == null) {
                bonds = mol.getConnectedBondsList(prev);
                cache.put(prev, bonds);
            }

            IBond bond = findBond(bonds, next, prev);
            if (bond == null)
                throw new IllegalStateException("FATAL - Atoms in patch were connected?");
            buffer.append(getBondSymbol(bond));
            buffer.append(getAtomSymbol(next));
            prev = next;
        }
        return buffer.toString();
    }

    private String encodePath(List<IAtom> apath, List<IBond> bpath, StringBuilder buffer) {
        buffer.setLength(0);
        IAtom prev = apath.get(0);
        buffer.append(getAtomSymbol(prev));
        for (int i = 1; i < apath.size(); i++) {
            final IAtom next  = apath.get(i);
            final IBond bond  = bpath.get(i-1);
            buffer.append(getBondSymbol(bond));
            buffer.append(getAtomSymbol(next));
        }
        return buffer.toString();
    }

    private int appendHash(int hash, String str) {
        int len = str.length();
        for (int i = 0; i < len; i++)
            hash = 31 * hash + str.charAt(0);
        return hash;
    }

    private int hashPath(List<IAtom> apath, List<IBond> bpath) {
        int hash = 0;
        hash = appendHash(hash, getAtomSymbol(apath.get(0)));
        for (int i = 1; i < apath.size(); i++) {
            final IAtom next  = apath.get(i);
            final IBond bond  = bpath.get(i-1);
            hash = appendHash(hash, getBondSymbol(bond));
            hash = appendHash(hash, getAtomSymbol(next));
        }
        return hash;
    }

    private int hashRevPath(List<IAtom> apath, List<IBond> bpath) {
        int hash = 0;
        int last = apath.size() - 1;
        hash = appendHash(hash, getAtomSymbol(apath.get(last)));
        for (int i = last-1; i >= 0; i--) {
            final IAtom next  = apath.get(i);
            final IBond bond  = bpath.get(i);
            hash = appendHash(hash, getBondSymbol(bond));
            hash = appendHash(hash, getAtomSymbol(next));
        }
        return hash;
    }

    private static final class State {
        private int    numPaths = 0;
        private Random rand = new Random();
        //private BitSet fp;
        private int[] fp;
        private IAtomContainer mol;
        private Set<IAtom> visited = new HashSet<>();
        private List<IAtom> apath = new ArrayList<>();
        private List<IBond> bpath = new ArrayList<>();
        private final int maxDepth;
        private final int fpsize;
        private Map<IAtom,List<IBond>> cache = new IdentityHashMap<>();
        public StringBuilder buffer = new StringBuilder();

        public State(IAtomContainer mol, int[] fp, int fpsize, int maxDepth) {
            this.mol = mol;
            this.fp  = fp;
            this.fpsize = fpsize;
            this.maxDepth = maxDepth;
        }

        List<IBond> getBonds(IAtom atom) {
            List<IBond> bonds = cache.get(atom);
            if (bonds == null) {
                bonds = mol.getConnectedBondsList(atom);
                cache.put(atom, bonds);
            }
            return bonds;
        }

        boolean visit(IAtom a) {
            return visited.add(a);
        }

        boolean unvisit(IAtom a) {
            return visited.remove(a);
        }

        void push(IAtom atom, IBond bond) {
            apath.add(atom);
            if (bond != null)
                bpath.add(bond);
        }

        void pop() {
            if (!apath.isEmpty())
                apath.remove(apath.size()-1);
            if (!bpath.isEmpty())
                bpath.remove(bpath.size()-1);
        }

        void addHash(int x) {
            numPaths++;
            rand.setSeed(x);
            // XXX: fp.set(x % size); would work just as well but would encode a
            //      different bit
            int idx = rand.nextInt(fpsize);
            fp[idx]++;
        }
    }

    private void traversePaths(State state, IAtom beg, IBond prev) throws CDKException {
        if (!hashPseudoAtoms && isPseudo(beg))
            return;
        state.push(beg, prev);
        state.addHash(encodeUniquePath(state.apath, state.bpath, state.buffer));
        if (state.numPaths > pathLimit)
            throw new CDKException("Too many paths! Structure is likely a cage, reduce path length or increase path limit");
        if (state.apath.size() < state.maxDepth) {
            for (IBond bond : state.getBonds(beg)) {
                if (bond.equals(prev))
                    continue;
                final IAtom nbr = bond.getOther(beg);
                if (state.visit(nbr)) {
                    traversePaths(state, nbr, bond);
                    state.unvisit(nbr); // traverse all paths
                }
            }
        }
        state.pop();
    }

    /**
     * Get all paths of lengths 0 to the specified length.
     *
     * This method will find all paths up to length N starting from each
     * atom in the molecule and return the unique set of such paths.
     *
     * @param container The molecule to search
     * @param searchDepth The maximum path length desired
     * @return A Map of path strings, keyed on themselves
     * @deprecated Use {@link #encodePaths(IAtomContainer, int, BitSet, int)}
     */
    @Deprecated
    protected int[] findPathes(IAtomContainer container, int searchDepth) throws CDKException {

        Set<Integer> hashes = new HashSet<>();

        Map<IAtom, List<IBond>> cache = new HashMap<>();
        StringBuilder buffer = new StringBuilder();
        for (IAtom startAtom : container.atoms()) {
            List<List<IAtom>> p = PathTools.getLimitedPathsOfLengthUpto(container, startAtom, searchDepth, pathLimit);
            for (List<IAtom> path : p) {
                if (hashPseudoAtoms || !hasPseudoAtom(path))
                    hashes.add(encodeUniquePath(container, cache, path, buffer));
            }
        }

        int   pos = 0;
        int[] result = new int[hashes.size()];
        for (Integer hash : hashes)
            result[pos++] = hash;

        return result;
    }

    //protected void encodePaths(IAtomContainer mol, int depth, BitSet fp, int size) throws CDKException {
    protected void encodePaths(IAtomContainer mol, int depth, int[] fp, int size) throws CDKException {
        State state = new State(mol, fp, size, depth+1);
        for (IAtom atom : mol.atoms()) {
            state.numPaths = 0;
            state.visit(atom);
            traversePaths(state, atom, null);
            state.unvisit(atom);
        }
    }

    private static boolean isPseudo(IAtom a) {
        return getElem(a) == 0;
    }

    private static boolean hasPseudoAtom(Iterable<IAtom> path) {
        for (IAtom atom : path)
            if (isPseudo(atom))
                return true;
        return false;
    }

    private int encodeUniquePath(IAtomContainer container, Map<IAtom, List<IBond>> cache, List<IAtom> path, StringBuilder buffer) {
        if (path.size() == 1)
            return getAtomSymbol(path.get(0)).hashCode();
        String forward = encodePath(container, cache, path, buffer);
        Collections.reverse(path);
        String reverse = encodePath(container, cache, path, buffer);
        Collections.reverse(path);

        final int x;
        if (reverse.compareTo(forward) < 0)
            x = forward.hashCode();
        else
            x = reverse.hashCode();
        return x;
    }

    /**
     * Compares atom symbols lexicographical
     * @param a atom a
     * @param b atom b
     * @return comparison &lt;0 a is less than b, &gt;0 a is more than b
     */
    private int compare(IAtom a, IAtom b) {
        final int elemA = getElem(a);
        final int elemB = getElem(b);
        if (elemA == elemB)
            return 0;
        return getAtomSymbol(a).compareTo(getAtomSymbol(b));
    }

    /**
     * Compares bonds symbols lexicographical
     * @param a bond a
     * @param b bond b
     * @return comparison &lt;0 a is less than b, &gt;0 a is more than b
     */
    private int compare(IBond a, IBond b) {
        return getBondSymbol(a).compareTo(getBondSymbol(b));
    }

    /**
     * Compares a path of atoms with it's self to give the
     * lexicographically lowest traversal (forwards or backwards).
     * @param apath path of atoms
     * @param bpath path of bonds
     * @return &lt;0 forward is lower &gt;0 reverse is lower
     */
    private int compare(List<IAtom> apath, List<IBond> bpath) {
        int i    = 0;
        int len = apath.size();
        int j    = len - 1;
        int cmp = compare(apath.get(i), apath.get(j));
        if (cmp != 0)
            return cmp;
        i++;
        j--;
        while (j != 0) {
            cmp = compare(bpath.get(i-1), bpath.get(j));
            if (cmp != 0) return cmp;
            cmp = compare(apath.get(i), apath.get(j));
            if (cmp != 0) return cmp;
            i++;
            j--;
        }
        return 0;
    }

    private int encodeUniquePath(List<IAtom> apath, List<IBond> bpath, StringBuilder buffer) {
        if (bpath.size() == 0)
            return getAtomSymbol(apath.get(0)).hashCode();
        final int x;
        if (compare(apath, bpath) >= 0) {
            x = hashPath(apath, bpath);
        } else {
            x = hashRevPath(apath, bpath);
        }
        return x;
    }

    private static int getElem(IAtom atom) {
        Integer elem = atom.getAtomicNumber();
        if (elem == null)
            elem = 0;
        return elem;
    }

    private String getAtomSymbol(IAtom atom) {
        // XXX: backwards compatibility
        // This is completely random, I believe the intention is because
        // paths were reversed with string manipulation to de-duplicate
        // (only the lowest lexicographically is stored) however this
        // doesn't work with multiple atom symbols:
        // e.g. Fe-C => C-eF vs C-Fe => eF-C
        // A dirty hack is to replace "common" symbols with single letter
        // equivalents so the reversing is less wrong
        switch (getElem(atom)) {
            case 0:  // *
                return "*";
            case 6:  // C
                return "C";
            case 7:  // N
                return "N";
            case 8:  // O
                return "O";
            case 17: // Cl
                return "X";
            case 35: // Br
                return "Z";
            case 14: // Si
                return "Y";
            case 33: // As
                return "D";
            case 3: // Li
                return "L";
            case 34: // Se
                return "E";
            case 11:  // Na
                return "G";
            case 20:  // Ca
                return "J";
            case 13:  // Al
                return "A";
        }
        return atom.getSymbol();
    }

    /**
     *  Gets the bondSymbol attribute of the Fingerprinter class
     *
     *@param  bond  Description of the Parameter
     *@return       The bondSymbol value
     */
    protected String getBondSymbol(IBond bond) {
        if (bond.isAromatic())
            return ":";
        switch (bond.getOrder()) {
            case SINGLE:
                return "-";
            case DOUBLE:
                return "=";
            case TRIPLE:
                return "#";
            default:
                return "";
        }
    }

    public void setPathLimit(int limit) {
        this.pathLimit = limit;
    }

    public void setHashPseudoAtoms(boolean value) {
        this.hashPseudoAtoms = value;
    }

    public int getSearchDepth() {
        return searchDepth;
    }

    @Override
    public int getSize() {
        return size;
    }

    @Override
    public ICountFingerprint getCountFingerprint(IAtomContainer container) throws CDKException {
        throw new UnsupportedOperationException();
    }

}
