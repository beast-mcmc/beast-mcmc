/*
 * SitePatterns.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.evolution.alignment;

import dr.evolution.datatype.AminoAcids;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.util.Pair;

import java.util.*;

import static dr.evolution.alignment.SitePatterns.CompressionType.*;

/**
 * Stores a set of site patterns. This differs from the simple Patterns
 * class because it stores the pattern index for each site. Thus it has
 * a connection to a single alignment.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class SitePatterns implements SiteList, dr.util.XHTMLable {
    private final boolean DEBUG = false;

    public enum CompressionType {
        // no compression - all patterns in order
        UNCOMPRESSED,
        // compressed - only unique patterns stored with weights
        UNIQUE_ONLY,
        // compressed - unique patterns and ambiguous matches stored with weights
        AMBIGUOUS_UNIQUE,
        // compressed - unique patterns and ambiguous constant patterns stored with weights
        AMBIGUOUS_CONSTANT
    }

    public static final CompressionType DEFAULT_COMPRESSION_TYPE = UNIQUE_ONLY;

    public static final double DEFAULT_AMBIGUITY_THRESHOLD = 0.25;

    public static final int MINIMUM_UNAMBIGUOUS = 2;

    private final boolean isCompressed;

    /**
     * the source alignment
     */
    protected final SiteList siteList;

    /**
     * number of sites
     */
    protected int siteCount = 0;

    /**
     * number of patterns
     */
    protected int patternCount = 0;

    /**
     * site -> site pattern
     */
    protected int[] sitePatternIndices;

    /**
     * count of invariant patterns
     */
    protected int invariantCount;

    /**
     * weights of each site pattern
     */
    protected double[] weights;

    /**
     * site patterns [site pattern][taxon]
     */
    protected int[][] patterns;

    protected double[][][] uncertainPatterns;

    private boolean uncertainSites = false;

    /**
     * Constructor
     */
    public SitePatterns(Alignment alignment) {
        this(alignment, null, -1, -1, 1);
    }

    /**
     * Constructor
     */
    public SitePatterns(Alignment alignment, TaxonList taxa) {
        this(alignment, taxa, -1, -1, 1, true);
    }

    public SitePatterns(Alignment alignment, TaxonList taxa, CompressionType compressionType) {
        this(alignment, taxa, -1, -1, 1, true, compressionType);
    }

    /**
     * Constructor
     */
    public SitePatterns(Alignment alignment, int from, int to, int every) {
        this(alignment, null, from, to, every);
    }

    /**
     * Constructor
     */
    public SitePatterns(Alignment alignment, TaxonList taxa, int from, int to, int every) {
        this(alignment, taxa, from, to, every,true);
    }

    public SitePatterns(Alignment alignment, TaxonList taxa, int from, int to, int every, boolean strip) {
        this(alignment, taxa, from, to, every, strip, DEFAULT_COMPRESSION_TYPE);
    }

    public SitePatterns(Alignment alignment, TaxonList taxa, int from, int to, int every, boolean strip, CompressionType compressionType) {
        this(alignment, taxa, from, to, every, strip, null, compressionType, DEFAULT_AMBIGUITY_THRESHOLD);
    }

    /**
     *
     * @param alignment The alignment
     * @param taxa The list of taxa - can be a subset of those in the alignment
     * @param from the first site to be included in the pattern list (zero indexed)
     * @param to the last site to be included in the pattern list (inclusive)
     * @param every skip every over every X sites
     * @param strip whether to strip completely ambiguous/gapped sites
     * @param constantSiteCounts a vector of counts of constant sites for each state (for where the alignment only includes variable sites)
     * @param compression Type of pattern/weight compression to use
     */
    public SitePatterns(Alignment alignment, TaxonList taxa, int from, int to, int every, boolean strip,
                        int[] constantSiteCounts, CompressionType compression, double ambiguityThreshold ) {

        this.siteList = alignment;
        isCompressed = compression != UNCOMPRESSED;

        if (taxa != null) {
            SimpleAlignment a = new SimpleAlignment();

            for (int i = 0; i < alignment.getSequenceCount(); i++) {
                if (taxa.getTaxonIndex(alignment.getTaxonId(i)) != -1) {
                    a.addSequence(alignment.getSequence(i));
                }
            }

            alignment = a;
        }

        if (constantSiteCounts != null) {
            if (constantSiteCounts.length != alignment.getStateCount()) {
                throw new IllegalArgumentException("Constant site count array length doesn't equal the number of states");
            }
        }

        addPatterns(alignment, from, to, every, strip, constantSiteCounts, compression, ambiguityThreshold);
    }

    /**
     * Constructor
     */
    public SitePatterns(SiteList siteList) {
        this(siteList, -1, -1, 1);
    }

    /**
     * Constructor
     */
    public SitePatterns(SiteList siteList, int from, int to, int every) {
        this(siteList, from, to, every, true, DEFAULT_COMPRESSION_TYPE, DEFAULT_AMBIGUITY_THRESHOLD);
    }

    /**
     * Constructor
     */
    public SitePatterns(SiteList siteList, int from, int to, int every, boolean strip) {
        this(siteList, from, to, every, strip, DEFAULT_COMPRESSION_TYPE, DEFAULT_AMBIGUITY_THRESHOLD);
    }

    /**
     * Constructor
     */
    public SitePatterns(SiteList siteList, int from, int to, int every, boolean strip, CompressionType compression, double ambiguityThreshold) {
        this.siteList = siteList;
        isCompressed = compression != UNCOMPRESSED;
        addPatterns(siteList, from, to, every, strip, null, compression, ambiguityThreshold);
    }

    public SiteList getSiteList() {
        return siteList;
    }

    /**
     * adds a set of patterns to the patternlist
     */
    private void addPatterns(SiteList siteList, int from, int to, int every, boolean strip, int[] constantSiteCounts,
                             CompressionType compression, double ambiguityThreshold) {
        if (siteList == null) {
            return;
        }

        if (from <= -1)
            from = 0;

        if (to <= -1)
            to = siteList.getSiteCount() - 1;

        if (every <= 0)
            every = 1;

        siteCount = ((to - from) / every) + 1;

        if (compression == AMBIGUOUS_CONSTANT || compression == AMBIGUOUS_UNIQUE) {
            // add some space for the set of constant sites
            siteCount += siteList.getStateCount();
        }

        patternCount = 0;

        patterns = new int[siteCount][];

        sitePatternIndices = new int[siteCount];

        weights = new double[siteCount];

        invariantCount = 0;

        uncertainSites = siteList.areUncertain();

        if (uncertainSites) {
            uncertainPatterns = new double[siteCount][][];
        }

        if (DEBUG) {
            System.err.println("Creating SitePatterns using compression type: " + compression.toString());
        }

        if (compression == AMBIGUOUS_CONSTANT || compression == AMBIGUOUS_UNIQUE) {
            // if the patterns are to be compressed then create the constant sites initially
            for (int i = 0; i < siteList.getStateCount(); i++) {
                int[] pattern = new int[siteList.getPatternLength()];
                for (int j = 0; j < siteList.getPatternLength(); j++) {
                    pattern[j] = i;
                }
                addPattern(pattern, constantSiteCounts != null ? constantSiteCounts[i] : 0, null);
            }
        }

        int site = 0;
        int count = 0;

        for (int i = from; i <= to; i += every) {
            int[] pattern = siteList.getSitePattern(i);
            double weight = siteList.getPatternWeight(i);

            if (uncertainSites) {
                sitePatternIndices[site] = addUncertainPattern(pattern, weight, siteList.getUncertainSitePattern(i));
            } else{
                // @todo - what is `strip` being used for?
                if (!strip || !isInvariant(pattern, false) ||
                        (!isGapped(pattern) &&
                                !isAmbiguous(pattern) &&
                                !isUnknown(pattern))) {

                    sitePatternIndices[site] = addPattern(pattern, weight, compression);

                    count += 1;
                } else {
                    sitePatternIndices[site] = -1;
                }
            }
            site++;
        }

        if (DEBUG) {
            System.err.println("Added " + count + " site patterns");

            if (compression != UNCOMPRESSED) {
                System.err.println("Compressed to " + patternCount + " unique patterns");
//                System.err.println("Constant patterns:");
//                for (int i = 0; i < getStateCount(); i++) {
//                    System.err.println("State " + getDataType().getCode(i) + ": " + weights[i]);
//                }

                System.err.println("All patterns:");
                for (int i = 0; i < getPatternCount(); i++) {
                    System.err.print("Pattern: ");
                    for (int j = 0; j < getPatternLength(); j++) {
                        System.err.print(getDataType().getCode(patterns[i][j]));
                    }
                    System.err.println(" - " + weights[i]);
                }
            }
        }


        if (compression != UNCOMPRESSED && compression != UNIQUE_ONLY) {
//            sortPatternsByWeight();
            compressAmbiguousPatterns(compression == AMBIGUOUS_CONSTANT, ambiguityThreshold);

            // these are no longer valid...
            sitePatternIndices = null;

            if (DEBUG) {
                System.err.println("Further compressed to " + patternCount + " ambiguously unique patterns");
//                System.err.println("Constant patterns:");
//                for (int i = 0; i < getStateCount(); i++) {
//                    System.err.println("State " + getDataType().getCode(i) + ": " + weights[i]);
//                }
                System.err.println("All patterns:");
                for (int i = 0; i < getPatternCount(); i++) {
                    System.err.print("Pattern: ");
                    for (int j = 0; j < getPatternLength(); j++) {
                        System.err.print(getDataType().getCode(patterns[i][j]));
                    }
                    System.err.println(" - " + weights[i]);
                }
            }
        }

        countInvariantSites();
    }

    /**
     * adds a pattern to the pattern list with the given weight
     *
     * @return the index of the pattern in the pattern list
     */
    private int addPattern(int[] pattern, double weight, CompressionType compression) {

        if (compression != UNCOMPRESSED) {
            // this will compress unique patterns, further compression of ambiguously similar
            // patterns is done in a later step
            for (int i = 0; i < patternCount; i++) {
                if (comparePatterns(patterns[i], pattern, false)) {
                    patterns[i] = pattern;
                    weights[i] += weight;
                    return i;
                }
            }
        }

        // new pattern - add it
        int index = patternCount;
        patterns[index] = pattern;
        weights[index] = weight;

        patternCount++;

        return index;
    }

    /**
     * adds an uncertain pattern (a matrix of probabilities) to the pattern list with the given weight
     *
     * @return the index of the pattern in the pattern list
     */
    private int addUncertainPattern(int[] pattern, double weight, double[][] uncertainty) {

        int index = patternCount;
        patterns[index] = pattern;
        weights[index] = weight;

        if (uncertainSites) {
            if (uncertainty == null) {
                uncertainPatterns[index] = new double[pattern.length][];
                for (int taxon = 0; taxon < pattern.length; ++taxon) {
                    double[] prob = new double[getDataType().getStateCount()];
                    int[] possibleStates = getDataType().getStates(pattern[taxon]);
                    for (int state : possibleStates) {
                        prob[state] = 1.0;
                    }
                    uncertainPatterns[index][taxon] = prob;
                }
            } else {
                uncertainPatterns[index] = uncertainty;
            }
        }

        patternCount++;

        return index;
    }

    private void countInvariantSites() {
        this.invariantCount = 0;
        for (int i = 0; i < patternCount; i++) {
            if (isInvariant(patterns[i], false)) {
                this.invariantCount += (int)weights[i];
            }
        }
    }

    /**
     * Orders the patterns by descending weight - for most data sets this will put the constant
     * patterns first.
     */
    private void sortPatternsByWeight() {
        List<Pair<int[], Double> > sortedPatterns = new ArrayList<>();
        for (int i = 0; i < patterns.length; i++) {
            sortedPatterns.add(new Pair<>(patterns[i], weights[i]));
        }
        // reverse sort (minus weight)
        sortedPatterns.sort(Comparator.comparingDouble(o -> -o.getSecond()));

        int i = 0;
        for (Pair<int[], Double> p : sortedPatterns) {
            patterns[i] = p.getFirst();
            weights[i] = p.getSecond();
            i++;
        }
    }

    private void compressAmbiguousPatterns(boolean constantOnly, double ambiguityThreshold) {
        int minimumUnambiguous = (int)((1.0 - ambiguityThreshold) * getPatternLength());
        minimumUnambiguous = Math.min(minimumUnambiguous, 2);

        // the first stateCount patterns are the constant ones
        for (int i = getStateCount(); i < patternCount; i++) {
            int count = constantOnly ? getStateCount() : i;
            for (int j = 0; j < count; j++) {
                if (patterns[j] != null) {
                    // the pattern should have at least 2 non-ambiguous characters
                    if (getCanonicalStateCount(patterns[i]) >= minimumUnambiguous &&
                            comparePatterns(patterns[i], patterns[j], true)) {
                        if (!constantOnly && getCanonicalStateCount(patterns[i]) > getCanonicalStateCount(patterns[j])) {
                            // if this is a less ambiguous pattern then this becomes the 'type' pattern
                            patterns[j] = patterns[i];
                        }
                        weights[j] += weights[i];
                        weights[i] = 0.0;
                        patterns[i] = null;
                        break;
                    }
                }
            }
        }

        // compress out all the nulls
        int i = 0;
        while (i < patternCount) {
            if (patterns[i] == null) {
                int j = i + 1;
                while (j < patternCount && patterns[j] == null) {
                    j += 1;
                }
                if (j < patternCount) {
                    patterns[i] = patterns[j];
                    weights[i] = weights[j];
                    patterns[j] = null;
                    weights[j] = 0.0;
                } else {
                    patternCount = i;
                }
            }
            i++;
        }
    }

    /**
     * @return true if the pattern contains a gap state
     */
    private boolean isGapped(int[] pattern) {
        int len = pattern.length;

        for (int i = 0; i < len; i++) {
            if (getDataType().isGapState(pattern[i])) {
                return true;
            }
        }
        return false;
    }

    /**
     * @return true if the pattern contains an ambiguous state
     */
    private boolean isAmbiguous(int[] pattern) {
        int len = pattern.length;

        for (int i = 0; i < len; i++) {
            if (getDataType().isAmbiguousState(pattern[i])) {
                return true;
            }
        }
        return false;
    }

    /**
     * @return true if the pattern contains an unknown state
     */
    private boolean isUnknown(int[] pattern) {
        int len = pattern.length;

        for (int i = 0; i < len; i++) {
            if (getDataType().isUnknownState(pattern[i])) {
                return true;
            }
        }
        return false;
    }

    /**
     * @return true if the pattern is invariant
     */
    private boolean isInvariant(int[] pattern, boolean ignoreAmbiguity) {
        int firstState = pattern[0];
        for (int i = 1; i < pattern.length; i++) {
            if (ignoreAmbiguity) {
                if (getDataType().areUnambiguouslyDifferent(firstState, pattern[i])) {
                    return false;
                }
            } else if (firstState != pattern[i]) {
                return false;
            }
        }

        return true;
    }

    /**
     * returns the number of canonical (non-ambiguous) states in the pattern
     * @param pattern
     * @return
     */
    private int getCanonicalStateCount(int[] pattern) {
        int count = 0;
        for (int state : pattern) {
            if (state < getDataType().getStateCount()) {
                count += 1;
            }
        }
        return count;
    }

    /**
     * compares two patterns
     *
     * @return true if they are identical
     */
    protected boolean comparePatterns(int[] pattern1, int[] pattern2) {
        return comparePatterns(pattern1, pattern2, false);
    }

    /**
     * compares two patterns
     *
     * @return true if they are identical
     */
    protected boolean comparePatterns(int[] pattern1, int[] pattern2, boolean allowAmbiguities) {

        if (!allowAmbiguities) {
            return Arrays.equals(pattern1, pattern2);

//            for (int i = 0; i < len; i++) {
//                if (pattern1[i] != pattern2[i]) {
//                    return false;
//                }
//            }
        } else {
            for (int i = 0; i < pattern1.length; i++) {
                if (getDataType().areUnambiguouslyDifferent(pattern1[i], pattern2[i])) {
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * @return number of invariant sites (these will be first in the list).
     */
    public int getInvariantCount() {
        return invariantCount;
    }

    // **************************************************************
    // SiteList IMPLEMENTATION
    // **************************************************************

    /**
     * @return number of sites
     */
    public int getSiteCount() {
        return siteCount;
    }

    /**
     * Gets the pattern of site as an array of state numbers (one per sequence)
     *
     * @return the site pattern at siteIndex
     */
    public int[] getSitePattern(int siteIndex) {
        final int sitePatternIndice = sitePatternIndices[siteIndex];
        return sitePatternIndice >= 0 ? patterns[sitePatternIndice] : null;
    }

    @Override
    public double[][] getUncertainSitePattern(int siteIndex) {
        throw new UnsupportedOperationException("getUncertainSitePattern not implemented yet");
    }

    /**
     * Gets the pattern index at a particular site
     *
     * @return the patternIndex
     */
    public int getPatternIndex(int siteIndex) {
        return sitePatternIndices[siteIndex];
    }

    /**
     * @return the sequence state at (taxon, site)
     */
    public int getState(int taxonIndex, int siteIndex) {
        final int sitePatternIndice = sitePatternIndices[siteIndex];
        // is that right?
        return sitePatternIndice >= 0 ? patterns[sitePatternIndice][taxonIndex] : getDataType().getGapState();
    }

    @Override
    public double[] getUncertainState(int taxonIndex, int siteIndex) {
        throw new UnsupportedOperationException("getUncertainState not implemented yet");
    }

    // **************************************************************
    // PatternList IMPLEMENTATION
    // **************************************************************

    /**
     * @return number of patterns
     */
    public int getPatternCount() {
        return patternCount;
    }

    /**
     * @return number of states for this siteList
     */
    public int getStateCount() {
        if (siteList == null) throw new RuntimeException("SitePatterns has no alignment");
        return siteList.getStateCount();
    }

    /**
     * Gets the length of the pattern strings which will usually be the
     * same as the number of taxa
     *
     * @return the length of patterns
     */
    public int getPatternLength() {
        return getTaxonCount();
    }

    /**
     * Gets the pattern as an array of state numbers (one per sequence)
     *
     * @return the pattern at patternIndex
     */
    public int[] getPattern(int patternIndex) {
        return patterns[patternIndex];
    }

    @Override
    public double[][] getUncertainPattern(int patternIndex) {
        throw new UnsupportedOperationException("getUncertainPattern not implemented yet");
    }

    /**
     * @return state at (taxonIndex, patternIndex)
     */
    public int getPatternState(int taxonIndex, int patternIndex) {
        return patterns[patternIndex][taxonIndex];
    }

    @Override
    public double[] getUncertainPatternState(int taxonIndex, int patternIndex) {
        return uncertainPatterns[patternIndex][taxonIndex];
    }

    /**
     * Gets the weight of a site pattern
     */
    public double getPatternWeight(int patternIndex) {
        return weights[patternIndex];
    }

    /**
     * @return the array of pattern weights
     */
    public double[] getPatternWeights() {
        return weights;
    }

    /**
     * @return the DataType of this siteList
     */
    public DataType getDataType() {
        if (siteList == null) throw new RuntimeException("SitePatterns has no alignment");
        return siteList.getDataType();
    }

    /**
     * @return the frequency of each state
     */
    public double[] getStateFrequencies() {
        return PatternList.Utils.empiricalStateFrequencies(this);
    }

    public boolean areUnique() {
        return isCompressed;
    }

    @Override
    public boolean areUncertain() {
        return uncertainSites;
    }

    // **************************************************************
    // TaxonList IMPLEMENTATION
    // **************************************************************

    /**
     * @return a count of the number of taxa in the list.
     */
    public int getTaxonCount() {
        if (siteList == null) throw new RuntimeException("SitePatterns has no alignment");
        return siteList.getTaxonCount();
    }

    /**
     * @return the ith taxon.
     */
    public Taxon getTaxon(int taxonIndex) {
        if (siteList == null) throw new RuntimeException("SitePatterns has no alignment");
        return siteList.getTaxon(taxonIndex);
    }

    /**
     * @return the ID of the ith taxon.
     */
    public String getTaxonId(int taxonIndex) {
        if (siteList == null) throw new RuntimeException("SitePatterns has no alignment");
        return siteList.getTaxonId(taxonIndex);
    }

    /**
     * returns the index of the taxon with the given id.
     */
    public int getTaxonIndex(String id) {
        if (siteList == null) throw new RuntimeException("SitePatterns has no alignment");
        return siteList.getTaxonIndex(id);
    }

    /**
     * returns the index of the given taxon.
     */
    public int getTaxonIndex(Taxon taxon) {
        if (siteList == null) throw new RuntimeException("SitePatterns has no alignment");
        return siteList.getTaxonIndex(taxon);
    }

    public List<Taxon> asList() {
        List<Taxon> taxa = new ArrayList<Taxon>();
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            taxa.add(getTaxon(i));
        }
        return taxa;
    }

    public Iterator<Taxon> iterator() {
        return new Iterator<Taxon>() {
            private int index = -1;

            public boolean hasNext() {
                return index < getTaxonCount() - 1;
            }

            public Taxon next() {
                index ++;
                return getTaxon(index);
            }

            public void remove() { /* do nothing */ }
        };
    }

    /**
     * @param taxonIndex the index of the taxon whose attribute is being fetched.
     * @param name       the name of the attribute of interest.
     * @return an object representing the named attributed for the given taxon.
     */
    public Object getTaxonAttribute(int taxonIndex, String name) {
        if (siteList == null) throw new RuntimeException("SitePatterns has no alignment");
        return siteList.getTaxonAttribute(taxonIndex, name);
    }

    // **************************************************************
    // Identifiable IMPLEMENTATION
    // **************************************************************

    protected String id = null;

    /**
     * @return the id.
     */
    public String getId() {
        return id;
    }

    /**
     * Sets the id.
     */
    public void setId(String id) {
        this.id = id;
    }

    // **************************************************************
    // XHTMLable IMPLEMENTATION
    // **************************************************************

    public String toXHTML() {
        String xhtml = "<p><em>Pattern List</em>  pattern count = ";
        xhtml += getPatternCount();
        xhtml += "  invariant count = ";
        xhtml += getInvariantCount();
        xhtml += "</p>";

        xhtml += "<pre>";

        int count, state;
        int type = getDataType().getType();

        count = getPatternCount();

        int length, maxLength = 0;
        for (int i = 0; i < count; i++) {
            length = Integer.toString((int) getPatternWeight(i)).length();
            if (length > maxLength)
                maxLength = length;
        }

        for (int i = 0; i < count; i++) {
            length = Integer.toString(i + 1).length();
            for (int j = length; j < maxLength; j++)
                xhtml += " ";
            xhtml += Integer.toString(i + 1) + ": ";

            length = Integer.toString((int) getPatternWeight(i)).length();
            xhtml += Integer.toString((int) getPatternWeight(i));
            for (int j = length; j <= maxLength; j++)
                xhtml += " ";

            for (int j = 0; j < getTaxonCount(); j++) {
                state = getPatternState(j, i);

                if (type == DataType.NUCLEOTIDES) {
                    xhtml += Nucleotides.INSTANCE.getChar(state) + " ";
                } else if (type == DataType.CODONS) {
                    xhtml += Codons.UNIVERSAL.getTriplet(state) + " ";
                } else {
                    xhtml += AminoAcids.INSTANCE.getChar(state) + " ";
                }
            }
            xhtml += "\n";
        }
        xhtml += "</pre>";
        return xhtml;
    }
}
