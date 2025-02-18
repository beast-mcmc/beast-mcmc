/*
 * BeastGenerator.java
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

package dr.app.beauti.generator;

import dr.app.beast.BeastVersion;
import dr.app.beauti.BeautiApp;
import dr.app.beauti.BeautiFrame;
import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.*;
import dr.app.beauti.types.*;
import dr.app.beauti.util.XMLWriter;
import dr.app.util.Arguments;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evoxml.AlignmentParser;
import dr.evoxml.DateParser;
import dr.evoxml.TaxaParser;
import dr.evoxml.TaxonParser;
import dr.inferencexml.model.CompoundLikelihoodParser;
import dr.inferencexml.operators.SimpleOperatorScheduleParser;
import dr.util.Attribute;
import dr.util.Pair;
import dr.util.Version;
import dr.xml.AttributeParser;
import dr.xml.XMLParser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * This class holds all the data for the current BEAUti Document
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class BeastGenerator extends Generator {

    private final static Version VERSION = new BeastVersion();
    private static final String MESSAGE_CAL_YULE = "Calibrated Yule requires 1 calibrated internal node \n" +
            "with a proper prior and monophyly enforced for each tree.";
    private final String MESSAGE_CAL = "\nas another element (taxon, sequence, taxon set, species, etc.):\nAll ids should be unique.";

    private final AlignmentGenerator alignmentGenerator;
    private final PatternListGenerator patternListGenerator;
    private final TreePriorGenerator treePriorGenerator;
    private final TreeLikelihoodGenerator treeLikelihoodGenerator;
    private final SubstitutionModelGenerator substitutionModelGenerator;
    private final InitialTreeGenerator initialTreeGenerator;
    private final TreeModelGenerator treeModelGenerator;
    private final ClockModelGenerator clockModelGenerator;
    private final OperatorsGenerator operatorsGenerator;
    private final ParameterPriorGenerator parameterPriorGenerator;
    private final LogGenerator logGenerator;
    private final TMRCAStatisticsGenerator tmrcaStatisticsGenerator;

    public BeastGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);

        alignmentGenerator = new AlignmentGenerator(options, components);
        patternListGenerator = new PatternListGenerator(options, components);
        tmrcaStatisticsGenerator = new TMRCAStatisticsGenerator(options, components);
        substitutionModelGenerator = new SubstitutionModelGenerator(options, components);
        treePriorGenerator = new TreePriorGenerator(options, components);
        treeLikelihoodGenerator = new TreeLikelihoodGenerator(options, components);

        initialTreeGenerator = new InitialTreeGenerator(options, components);
        treeModelGenerator = new TreeModelGenerator(options, components);
        clockModelGenerator = new ClockModelGenerator(options, components);

        operatorsGenerator = new OperatorsGenerator(options, components);
        parameterPriorGenerator = new ParameterPriorGenerator(options, components);
        logGenerator = new LogGenerator(options, components);

    }

    /**
     * Checks various options to check they are valid. Throws IllegalArgumentExceptions with
     * descriptions of the problems.
     *
     * @throws IllegalArgumentException if there is a problem with the current settings
     */
    public void checkOptions() throws GeneratorException {
        try{

            //++++++++++++++++ Taxon List ++++++++++++++++++
            TaxonList taxonList = options.taxonList;
            Set<String> ids = new HashSet<String>();

            ids.add(TaxaParser.TAXA);
            ids.add(AlignmentParser.ALIGNMENT);

            if (taxonList != null) {
                if (taxonList.getTaxonCount() < 2) {
                    throw new GeneratorException("BEAST requires at least two taxa to run.");
                }

                for (int i = 0; i < taxonList.getTaxonCount(); i++) {
                    Taxon taxon = taxonList.getTaxon(i);
                    if (ids.contains(taxon.getId())) {
                        throw new GeneratorException("A taxon has the same id," + taxon.getId() +  MESSAGE_CAL);
                    }
                    ids.add(taxon.getId());
                }
            }

            //++++++++++++++++ Taxon Sets ++++++++++++++++++
            for (PartitionTreeModel model : options.getPartitionTreeModels()) {
                // should be only 1 calibrated internal node with a proper prior and monophyletic for each tree at moment
                if (model.getPartitionTreePrior().getNodeHeightPrior() == TreePriorType.YULE_CALIBRATION) {
                    if (options.treeModelOptions.isNodeCalibrated(model) < 0) // invalid node calibration
                        throw new GeneratorException(MESSAGE_CAL_YULE);

                    if (options.treeModelOptions.isNodeCalibrated(model) > 0) { // internal node calibration
                        List taxonSetsList = options.getKeysFromValue(options.taxonSetsTreeModel, model);
                        if (taxonSetsList.size() != 1 || !options.taxonSetsMono.get(taxonSetsList.get(0))) { // 1 tmrca per tree && monophyletic
                            throw new GeneratorException(MESSAGE_CAL_YULE, BeautiFrame.TAXON_SETS);
                        }
                    }
                }
            }

            for (Taxa taxa : options.taxonSets) {
                // AR - we should allow single taxon taxon sets...
                if (taxa.getTaxonCount() < 1 // && !options.taxonSetsIncludeStem.get(taxa)
                        ) {
                    throw new GeneratorException(
                            "Taxon set, " + taxa.getId() + ", should contain \n" +
                                    "at least one taxa. Please go back to Taxon Sets \n" +
                                    "panel to correct this.", BeautiFrame.TAXON_SETS);
                }
                if (ids.contains(taxa.getId())) {
                    throw new GeneratorException("A taxon set has the same id," + taxa.getId() +
                            MESSAGE_CAL, BeautiFrame.TAXON_SETS);
                }
                ids.add(taxa.getId());
            }

            //++++++++++++++++ Traits ++++++++++++++++++
            // missing data is not necessarily an issue...
//        for (TraitData trait : options.traits) {
//            for (int i = 0; i < trait.getTaxaCount(); i++) {
////                System.out.println("Taxon " + trait.getTaxon(i).getId() + " : [" + trait.getTaxon(i).getAttribute(trait.getName()) + "]");
//                if (!trait.hasValue(i))
//                    throw new IllegalArgumentException("Taxon " + trait.getTaxon(i).getId() +
//                            " has no value for Trait " + trait.getName());
//            }
//        }

            //++++++++++++++++ Tree Prior ++++++++++++++++++
//        if (options.isShareSameTreePrior()) {
            if (options.getPartitionTreeModels().size() > 1) { //TODO not allowed multi-prior yet
                for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
                    if (prior.getNodeHeightPrior() == TreePriorType.GMRF_SKYRIDE) {
                        throw new GeneratorException("For the Skyride, tree model/tree prior combination not implemented by BEAST." +
                                "\nThe Skyride is only available for a single tree model partition in this release.", BeautiFrame.TREES);
                    }
                }
            }

            for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
                if ((prior.getNodeHeightPrior() == TreePriorType.SKYGRID || prior.getNodeHeightPrior() == TreePriorType.SKYGRID_HMC)
                        && Double.isNaN(prior.getSkyGridInterval())) {
                    throw new GeneratorException("The Skygrid cut-off time must be set and greater than 0.0.", BeautiFrame.TREES);
                }
            }

            //+++++++++++++++ Starting tree ++++++++++++++++
            for (PartitionTreeModel model : options.getPartitionTreeModels()) {
                if (model.getStartingTreeType() == StartingTreeType.USER) {
                    if (model.getUserStartingTree() == null) {
                        throw new GeneratorException("Please select a starting tree in " + BeautiFrame.TREES + " panel, " +
                                "\nwhen choosing user specified starting tree option.", BeautiFrame.TREES);
                    }
                }
            }

            //++++++++++++++++ Random local clock model validation ++++++++++++++++++
            for (PartitionClockModel model : options.getPartitionClockModels()) {
                // 1 random local clock CANNOT have different tree models
                if (model.getClockType() == ClockType.RANDOM_LOCAL_CLOCK) { // || AUTOCORRELATED_LOGNORMAL
                    PartitionTreeModel treeModel = null;
                    for (AbstractPartitionData pd : options.getDataPartitions(model)) { // only the PDs linked to this tree model
                        if (treeModel != null && treeModel != pd.getPartitionTreeModel()) {
                            throw new GeneratorException("A single random local clock cannot be applied to multiple trees.", BeautiFrame.CLOCK_MODELS);
                        }
                        treeModel = pd.getPartitionTreeModel();
                    }
                }
            }

            //++++++++++++++++ Tree Model ++++++++++++++++++
            for (PartitionTreeModel model : options.getPartitionTreeModels()) {
                int numOfTaxa = -1;
                for (AbstractPartitionData pd : options.getDataPartitions(model)) {
                    if (pd.getTaxonCount() > 0) {
                        if (numOfTaxa > 0) {
                            if (numOfTaxa != pd.getTaxonCount()) {
                                throw new GeneratorException("Partitions with different taxa cannot share the same tree.", BeautiFrame.DATA_PARTITIONS);
                            }
                        } else {
                            numOfTaxa = pd.getTaxonCount();
                        }
                    }
                }
            }

            //+++++++++++++++++ ThorneyBEAST / clock model Verification ++++++++++++++

            for (PartitionClockModel model : options.getPartitionClockModels()) {
                for (AbstractPartitionData pd : options.getDataPartitions(model)) {
                    if(pd.getDataType().getType() ==DataType.TREE){
                        PartitionTreeModel tm = pd.getPartitionTreeModel();
                        if(tm.isUsingThorneyBEAST() && !PartitionTreeModel.thorneyCompatibleClocks.contains(model.getClockType())){
                            throw new GeneratorException("Gradients for "+model.getClockType()+ 
                            " are not implemented for Thorney BEAST." +
                            "Please use a clock model that does not rely on gradients.", BeautiFrame.CLOCK_MODELS);

                        }
                    }
                }
            }


            //++++++++++++++++ Prior Bounds ++++++++++++++++++
            for (Parameter param : options.selectParameters()) {
                if (!Double.isNaN(param.getInitial()) ) {
                    if (param.isTruncated && (param.getInitial() < param.truncationLower || param.getInitial() > param.truncationUpper)) {
                        throw new GeneratorException("Parameter \"" + param.getName() + "\":" +
                                "\ninitial value " + param.getInitial() + " is NOT in the range [" + param.truncationLower + ", " + param.truncationUpper + "]," +
                                "\nor this range is wrong. Please check the Prior panel.", BeautiFrame.PRIORS);
                    } else if (param.priorType == PriorType.UNIFORM_PRIOR && (param.getInitial() < param.uniformLower || param.getInitial() > param.uniformUpper)) {
                        throw new GeneratorException("Parameter \"" + param.getName() + "\":" +
                                "\ninitial value " + param.getInitial() + " is NOT in the range [" + param.uniformLower + ", " + param.uniformUpper + "]," +
                                "\nor this range is wrong. Please check the Prior panel.", BeautiFrame.PRIORS);
                    }
                    if (param.isNonNegative && param.getInitial() < 0.0) {
                        throw new GeneratorException("Parameter \"" + param.getName() + "\":" +
                                "\ninitial value " + param.getInitial() + " should be non-negative. Please check the Prior panel.", BeautiFrame.PRIORS);
                    }

                    if (param.isZeroOne && (param.getInitial() < 0.0 || param.getInitial() > 1.0)) {
                        throw new GeneratorException("Parameter \"" + param.getName() + "\":" +
                                "\ninitial value " + param.getInitial() + " should lie in the interval [0, 1]. Please check the Prior panel.", BeautiFrame.PRIORS);
                    }
                }
            }

            checkComponentOptions();

            // add other tests and warnings here
            // Speciation model with dated tips
            // Sampling rates without dated tips or priors on rate or nodes

        } catch (Exception e) {
            // catch any other exceptions here and rethrow to generate messages
            throw new GeneratorException(e.getMessage());
        }


    }
    public boolean checkUserTreeIsBifurcating(){
        for (PartitionTreeModel model : options.getPartitionTreeModels()) {
            if(model.getStartingTreeType()==StartingTreeType.USER){
                if(!(isBifurcatingTree(model.getUserStartingTree(),model.getUserStartingTree().getRoot()))){
                    return false;
                }
            }
        }
        return true;
    }
    public boolean isBifurcatingTree(Tree tree, NodeRef node) {
        if (tree.getChildCount(node) > 2) return false;
        for (int i = 0; i < tree.getChildCount(node); i++) {
            if (!isBifurcatingTree(tree, tree.getChild(node, i))) return false;
        }
        return true;
    }

    /**
     * Generate a beast xml file from these beast options
     *
     * @param file File
     * @throws java.io.IOException IOException
     * @throws dr.app.util.Arguments.ArgumentException
     *                             ArgumentException
     */
    public void generateXML(File file) throws GeneratorException, IOException, Arguments.ArgumentException {

        XMLWriter writer = new XMLWriter(new BufferedWriter(new FileWriter(file)));

        writer.writeText("<?xml version=\"1.0\" standalone=\"yes\"?>");
        writer.writeComment("Generated by BEAUTi " + VERSION.getVersionString(),
                "      by Alexei J. Drummond, Andrew Rambaut and Marc A. Suchard",
                "      Department of Computer Science, University of Auckland and",
                "      Institute of Evolutionary Biology, University of Edinburgh",
                "      David Geffen School of Medicine, University of California, Los Angeles",
                "      http://beast.community/");
        writer.writeOpenTag("beast", new Attribute.Default<String>("version", BeautiApp.VERSION.getVersion()) );
        writer.writeText("");

        // this gives any added implementations of the 'Component' interface a
        // chance to generate XML at this point in the BEAST file.
        generateInsertionPoint(ComponentGenerator.InsertionPoint.BEFORE_TAXA, writer);

        if (options.originDate != null) {
            // Create a dummy taxon whose job is to specify the origin date
            Taxon originTaxon = new Taxon("originTaxon");
            options.originDate.setUnits(options.units);
            originTaxon.setDate(options.originDate);
            writeTaxon(originTaxon, true, false, writer);
        }


        //++++++++++++++++ Taxon List ++++++++++++++++++
        try {
            // write complete taxon list
            writeTaxa(options.taxonList, writer);

            writer.writeText("");

        } catch (Exception e) {
            e.printStackTrace(System.err);
            throw new GeneratorException("Taxon list generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ Taxon Sets ++++++++++++++++++
        List<Taxa> taxonSets = options.taxonSets;
        try {
            if (taxonSets != null && !taxonSets.isEmpty()) {
                tmrcaStatisticsGenerator.writeTaxonSets(writer, taxonSets);
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Taxon sets generation has failed:\n" + e.getMessage());
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TAXA, writer);

        //++++++++++++++++ Alignments ++++++++++++++++++
        List<Alignment> alignments = new ArrayList<Alignment>();
        try {
            for (AbstractPartitionData partition : options.dataPartitions) {
                Alignment alignment = null;
                if (partition instanceof PartitionData) {
                    alignment = ((PartitionData) partition).getAlignment();
                }
                if (alignment != null && !alignments.contains(alignment)) {
                    alignments.add(alignment);
                }
            }
            if (alignments.size() > 0) {
                alignmentGenerator.writeAlignments(alignments, writer);
                generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_SEQUENCES, writer);
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Alignments generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ Pattern Lists ++++++++++++++++++
        try {
            // Construct pattern lists even if sampling from a null alignment
            //if (!options.samplePriorOnly) {
                for (AbstractPartitionData partition : options.dataPartitions) { // Each PD has one TreeLikelihood
                    if (partition.getTaxonList() != null) {
                        switch (partition.getDataType().getType()) {
                            case DataType.NUCLEOTIDES:
                            case DataType.AMINO_ACIDS:
                            case DataType.CODONS:
                            case DataType.COVARION:
                            case DataType.TWO_STATES:
                                patternListGenerator.writePatternList((PartitionData) partition, writer);
                                break;

                            case DataType.GENERAL:
                            case DataType.CONTINUOUS:
                                // no patternlist for trait data - discrete (general) data type uses an
                                // attribute patterns which is generated next bit of this method.
                                break;

                            case DataType.TREE:
                            case DataType.DUMMY:
                                //Do nothing
                                break;

                            default:
                                throw new IllegalArgumentException("Unsupported data type");
                        }
                        writer.writeText("");
                    }
                }
            //}
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Pattern lists generation has failed:\n" + e.getMessage());
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_PATTERNS, writer);

        //++++++++++++++++ Tree Prior Model ++++++++++++++++++
        try {
            for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
                treePriorGenerator.writeTreePriorModel(prior, writer);
                writer.writeText("");
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Tree prior model generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ Starting Tree ++++++++++++++++++
        try {
            for (PartitionTreeModel model : options.getPartitionTreeModels()) {
                initialTreeGenerator.writeStartingTree(model, writer);
                writer.writeText("");
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Starting tree generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ Tree Model +++++++++++++++++++
        try {
            for (PartitionTreeModel model : options.getPartitionTreeModels()) {
                treeModelGenerator.writeTreeModel(model, writer);
                writer.writeText("");
            }

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TREE_MODEL, writer);
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Tree model generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ Statistics ++++++++++++++++++
        try {
            if (taxonSets != null && taxonSets.size() > 0) {
                tmrcaStatisticsGenerator.writeTMRCAStatistics(writer);
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("TMRCA statistics generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ Tree Prior Likelihood ++++++++++++++++++
        try {
            for (PartitionTreeModel model : options.getPartitionTreeModels()) {
                treePriorGenerator.writePriorLikelihood(model, writer);
                writer.writeText("");
            }

            for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
                treePriorGenerator.writeMultiLociTreePriors(prior, writer);
            }

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TREE_PRIOR, writer);
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Tree prior likelihood generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ Branch Rates Model ++++++++++++++++++
        try {
            for (PartitionClockModel model : options.getPartitionClockModels()) {
                clockModelGenerator.writeBranchRatesModel(model, writer);
                writer.writeText("");
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Branch rates model generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ Substitution Model & Site Model ++++++++++++++++++
        try {
            for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
                substitutionModelGenerator.writeSubstitutionSiteModel(model, writer);
                writer.writeText("");
            }

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_SUBSTITUTION_MODEL, writer);
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Substitution model or site model generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ AllMus parameter ++++++++++++++++++
        try {
            for (PartitionClockModel model : options.getPartitionClockModels()) {
                clockModelGenerator.writeAllMus(model, writer);
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Clock model generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ Site Model ++++++++++++++++++
//        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
//            substitutionModelGenerator.writeSiteModel(model, writer); // site model
//            substitutionModelGenerator.writeAllMus(model, writer); // allMus
//            writer.writeText("");
//        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_SITE_MODEL, writer);

        //++++++++++++++++ Tree Likelihood ++++++++++++++++++
        try {
            Map<Pair<Pair<PartitionTreeModel, PartitionClockModel>, DataType>, List<PartitionData>> partitionLists = new HashMap<Pair<Pair<PartitionTreeModel, PartitionClockModel>, DataType>, List<PartitionData>>();
            options.multiPartitionLists.clear();
            options.otherPartitions.clear();

            for (AbstractPartitionData partition : options.dataPartitions) {
                // generate tree likelihoods for alignment data partitions
                if (partition.getTaxonList() != null) {

                    if (treeLikelihoodGenerator.canUseMultiPartition(partition)) {
                        // all sequence partitions of the same type as the first into the list for use in a
                        // MultipartitionTreeDataLikelihood. Must also share the same tree, clock model and not be doing
                        // ancestral reconstruction or counting
                        Pair<Pair<PartitionTreeModel, PartitionClockModel>, DataType> key = new Pair(new Pair(partition.getPartitionTreeModel(),partition.getPartitionClockModel()), partition.getDataType());
                        List<PartitionData> partitions = partitionLists.get(key);

                        if (partitions == null) {
                            partitions = new ArrayList<PartitionData>();
                            options.multiPartitionLists.add(partitions);
                        }

                        partitions.add((PartitionData) partition);
                        partitionLists.put(key, partitions);

                    } else {
                        options.otherPartitions.add(partition);
                    }


                }
            }

            treeLikelihoodGenerator.writeAllTreeLikelihoods(writer);

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TREE_LIKELIHOOD, writer);
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Tree likelihood generation has failed:\n" + e.getMessage());
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TRAITS, writer);

        //++++++++++++++++ Operators ++++++++++++++++++
        try {
            generateInsertionPoint(ComponentGenerator.InsertionPoint.BEFORE_OPERATORS, writer);

            List<Operator> operators = options.selectOperators();
            operatorsGenerator.writeOperatorSchedule(operators, writer);
            writer.writeText("");

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_OPERATORS, writer);
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Operators generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ MCMC ++++++++++++++++++
        try {
            // XMLWriter writer, List<PartitionSubstitutionModel> models,
            writeMCMC(writer);
            writer.writeText("");

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_MCMC, writer);
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("MCMC or log generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++  ++++++++++++++++++
        try {
            writeTimerReport(writer);
            writer.writeText("");
            if (options.performTraceAnalysis) {
                writeTraceAnalysis(writer);
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("The last part of XML generation has failed:\n" + e.getMessage());
        }

        writer.writeCloseTag("beast");
        writer.flush();
        writer.close();
    }

    /**
     * Generate a taxa block from these beast options
     *
     * @param writer    the writer
     * @param taxonList the taxon list to write
     * @throws dr.app.util.Arguments.ArgumentException
     *          ArgumentException
     */
    private void writeTaxa(TaxonList taxonList, XMLWriter writer) throws Arguments.ArgumentException {
        // -1 (single taxa), 0 (1st gene of multi-taxa)

        writer.writeComment("The list of taxa to be analysed (can also include dates/ages).",
                "ntax=" + taxonList.getTaxonCount());
        writer.writeOpenTag(TaxaParser.TAXA, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, TaxaParser.TAXA)});

        boolean hasAttr = options.traits.size() > 0;

        boolean firstDate = true;
        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            Taxon taxon = taxonList.getTaxon(i);

            boolean hasDate = false;
            if (options.useTipDates) {
                hasDate = TaxonList.Utils.hasAttribute(taxonList, i, dr.evolution.util.Date.DATE);
            }

            if (hasDate) {
                dr.evolution.util.Date date = (dr.evolution.util.Date) taxon.getAttribute(dr.evolution.util.Date.DATE);

                if (firstDate) {
                    options.units = date.getUnits();
                    firstDate = false;
                } else {
                    if (options.units != date.getUnits()) {
                        System.err.println("Error: Units in dates do not match.");
                    }
                }
            }

            writeTaxon(taxon, hasDate, hasAttr, writer);

        }

        writer.writeCloseTag(TaxaParser.TAXA);
    }

    /**
     * Generate a taxa block from these beast options
     *
     * @param writer    the writer
     * @param taxon the taxon to write
     * @throws dr.app.util.Arguments.ArgumentException
     *          ArgumentException
     */
    private void writeTaxon(Taxon taxon, boolean hasDate, boolean hasAttr, XMLWriter writer) throws Arguments.ArgumentException {

        writer.writeTag(TaxonParser.TAXON, new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, taxon.getId())},
                !(hasDate || hasAttr)); // false if any of hasDate or hasAttr is true


        if (hasDate) {
            dr.evolution.util.Date date = (dr.evolution.util.Date) taxon.getAttribute(dr.evolution.util.Date.DATE);

            Attribute[] attributes;
            if (date.getUncertainty() > 0.0) {
                attributes = new Attribute[] {
                        new Attribute.Default<Double>(DateParser.VALUE, date.getTimeValue()),
                        new Attribute.Default<String>(DateParser.DIRECTION, date.isBackwards() ? DateParser.BACKWARDS : DateParser.FORWARDS),
                        new Attribute.Default<String>(DateParser.UNITS, Units.Utils.getDefaultUnitName(options.units)),
                        new Attribute.Default<Double>(DateParser.UNCERTAINTY, date.getUncertainty())
                };
            } else {
                attributes = new Attribute[] {
                        new Attribute.Default<Double>(DateParser.VALUE, date.getTimeValue()),
                        new Attribute.Default<String>(DateParser.DIRECTION, date.isBackwards() ? DateParser.BACKWARDS : DateParser.FORWARDS),
                        new Attribute.Default<String>(DateParser.UNITS, Units.Utils.getDefaultUnitName(options.units))
                        //new Attribute.Default("origin", date.getOrigin()+"")
                };
            }

            writer.writeTag(dr.evolution.util.Date.DATE, attributes, true);
        }

        for (TraitData trait : options.traits) {
            // there is no harm in allowing the species trait to be listed in the taxa
//                if (!trait.getName().equalsIgnoreCase(TraitData.TRAIT_SPECIES)) {
            writer.writeOpenTag(AttributeParser.ATTRIBUTE, new Attribute[]{
                    new Attribute.Default<String>(Attribute.NAME, trait.getName())});

            // denotes missing data using '?'
            writer.writeText(taxon.containsAttribute(trait.getName()) ? taxon.getAttribute(trait.getName()).toString() : "?");
            writer.writeCloseTag(AttributeParser.ATTRIBUTE);
//                }
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_TAXON, taxon, writer);

        if (hasDate || hasAttr) writer.writeCloseTag(TaxonParser.TAXON);

    }

    /**
     * Write the timer report block.
     *
     * @param writer the writer
     */
    public void writeTimerReport(XMLWriter writer) {
        writer.writeOpenTag("report");
        writer.writeOpenTag("property", new Attribute.Default<String>("name", "timer"));
        writer.writeIDref("mcmc", "mcmc");
        writer.writeCloseTag("property");
        writer.writeCloseTag("report");
    }

    /**
     * Write the trace analysis block.
     *
     * @param writer the writer
     */
    public void writeTraceAnalysis(XMLWriter writer) {
        writer.writeTag(
                "traceAnalysis",
                new Attribute[]{
                        new Attribute.Default<String>("fileName", options.logFileName)
                },
                true
        );
    }

    /**
     * Write the MCMC block.
     *
     * @param writer XMLWriter
     */
    public void writeMCMC(XMLWriter writer) {
        writer.writeComment("Define MCMC");

        List<Attribute> attributes = new ArrayList<Attribute>();
        attributes.add(new Attribute.Default<String>(XMLParser.ID, "mcmc"));
        attributes.add(new Attribute.Default<Integer>("chainLength", options.chainLength));
        attributes.add(new Attribute.Default<String>("autoOptimize", options.autoOptimize ? "true" : "false"));

        if (options.operatorAnalysis) {
            attributes.add(new Attribute.Default<String>("operatorAnalysis", options.operatorAnalysisFileName));
        }

        writer.writeOpenTag("mcmc", attributes);

        if (options.hasData()) {
            writer.writeOpenTag(CompoundLikelihoodParser.JOINT, new Attribute.Default<String>(XMLParser.ID, "joint"));
        }

        // write prior block
        writer.writeOpenTag(CompoundLikelihoodParser.PRIOR, new Attribute.Default<String>(XMLParser.ID, "prior"));

        parameterPriorGenerator.writeParameterPriors(writer);

        for (PartitionTreeModel model : options.getPartitionTreeModels()) {
            PartitionTreePrior prior = model.getPartitionTreePrior();
            treePriorGenerator.writePriorLikelihoodReference(prior, model, writer);
            writer.writeText("");
        }

        for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
            treePriorGenerator.writeMultiLociLikelihoodReference(prior, writer);
            writer.writeText("");
        }

        clockModelGenerator.writeClockLikelihoodReferences(writer);

        generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_MCMC_PRIOR, writer);

        writer.writeCloseTag(CompoundLikelihoodParser.PRIOR);

        if (options.hasData()) {
            // write likelihood block
            writer.writeOpenTag(CompoundLikelihoodParser.LIKELIHOOD, new Attribute.Default<String>(XMLParser.ID, "likelihood"));

            treeLikelihoodGenerator.writeTreeLikelihoodReferences(writer);

            generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_MCMC_LIKELIHOOD, writer);

            writer.writeCloseTag(CompoundLikelihoodParser.LIKELIHOOD);

            writer.writeCloseTag(CompoundLikelihoodParser.JOINT);
        }

        writer.writeIDref(SimpleOperatorScheduleParser.OPERATOR_SCHEDULE, "operators");

        // write log to screen
        logGenerator.writeLogToScreen(writer, clockModelGenerator, substitutionModelGenerator);

        // write log to file
        logGenerator.writeLogToFile(writer, treePriorGenerator, clockModelGenerator,
                substitutionModelGenerator, treeLikelihoodGenerator, tmrcaStatisticsGenerator);

        // write tree log to file
        logGenerator.writeTreeLogToFile(writer);

        //write current state of chain to checkpoint file
        logGenerator.writeCheckpointToFile(writer);

        writer.writeCloseTag("mcmc");
    }

}
