/*
 * TreeModelGenerator.java
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

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionTreeModel;
import dr.app.beauti.util.XMLWriter;
import dr.evomodel.bigfasttree.thorney.ConstrainedTreeModel;
import dr.evomodel.bigfasttree.thorney.RootHeightProxyParameter;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.EmpiricalTreeDistributionModel;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightProxyParameter;
import dr.evomodelxml.bigfasttree.thorney.ConstrainedTreeModelParser;
import dr.evomodelxml.coalescent.OldCoalescentSimulatorParser;
import dr.evomodelxml.tree.*;
import dr.evoxml.NewickParser;
import dr.evoxml.TaxaParser;
import dr.evoxml.UPGMATreeParser;
import dr.inference.model.ParameterParser;
import dr.inference.model.Statistic;
import dr.inference.model.StatisticParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class TreeModelGenerator extends Generator {

    public TreeModelGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

    void writeTreeModel(PartitionTreeModel model, XMLWriter writer) {
        if (model.isUsingEmpiricalTrees()) {
                writeEmpiricalTreeModel(model, writer);
        } else if(model.isUsingThorneyBEAST()){
                writeConstrainedTreeModel(model, writer);
        }
        else {
            writeDefaultTreeModel(model, writer);
        }

    }

    /**
     * Write constrained tree model XML block.
     *
     * @param model
     * @param writer the writer
     */


    void writeConstrainedTreeModel(PartitionTreeModel model, XMLWriter writer) {
        // write newick with dates
        //simulate tree
        //write tree model

        String prefix = model.getPrefix();

        final String treeModelName = prefix + DefaultTreeModel.TREE_MODEL; // treemodel.treeModel or treeModel
        final String STARTING_TREE = InitialTreeGenerator.STARTING_TREE;
// newick written by starting tree 

        writer.writeComment("A constrained tree model");
        writer.writeTag(ConstrainedTreeModel.CONSTRAINED_TREE_MODEL,
                new Attribute[] {
                        new Attribute.Default<>(XMLParser.ID, treeModelName),
                }, false);
                writer.writeIDref(OldCoalescentSimulatorParser.COALESCENT_TREE, prefix + STARTING_TREE); 

                writer.writeTag(ConstrainedTreeModelParser.CONSTRAINTS_TREE,false);
                        writer.writeIDref(NewickParser.NEWICK, prefix + "constraintsTree" );// TODO magic values!
                writer.writeCloseTag(ConstrainedTreeModelParser.CONSTRAINTS_TREE);

                writer.writeComment("Parameter proxies for node heights");
                
                writer.writeOpenTag(TreeModelParser.ROOT_HEIGHT);
                        writer.writeTag(ParameterParser.PARAMETER,
                                new Attribute.Default<String>(XMLParser.ID, treeModelName + "." + TreeModelParser.ROOT_HEIGHT), true);
                writer.writeCloseTag(TreeModelParser.ROOT_HEIGHT);
                
                writer.writeOpenTag(TreeModelParser.NODE_HEIGHTS, new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"));
                        writer.writeTag(ParameterParser.PARAMETER,
                                new Attribute.Default<String>(XMLParser.ID, treeModelName + "." + "internalNodeHeights"), true);
                writer.writeCloseTag(TreeModelParser.NODE_HEIGHTS);
        
                writer.writeOpenTag(TreeModelParser.NODE_HEIGHTS,
                        new Attribute[]{
                                new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
                                new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "true")
                        });
                        writer.writeTag(ParameterParser.PARAMETER,
                                new Attribute.Default<String>(XMLParser.ID, treeModelName + "." + "allInternalNodeHeights"), true);
                writer.writeCloseTag(TreeModelParser.NODE_HEIGHTS);

        // Insertion point for leaf sampling (not currently implemented - parser will throw an error)
        generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_TREE_MODEL, model, writer);

        writer.writeCloseTag(ConstrainedTreeModel.CONSTRAINED_TREE_MODEL);


        writeTreeModelStatistics(model, writer);

    }

    /**
     * Write an empirical tree model XML block.
     *
     * @param model
     * @param writer the writer
    <empiricalTreeDistributionModel id="treeModel" fileName="subset1.trees">
    <taxa idref="subset1"/>
    </empiricalTreeDistributionModel>
    <statistic id="treeModel.currentTree" name="Current Tree">
    <empiricalTreeDistributionModel idref="treeModel"/>
    </statistic>
     */
    void writeEmpiricalTreeModel(PartitionTreeModel model, XMLWriter writer) {
        String prefix = model.getPrefix();

        final String treeModelName = prefix + DefaultTreeModel.TREE_MODEL; // treemodel.treeModel or treeModel

        writer.writeComment("An empirical distribution of trees");
        writer.writeTag(EmpiricalTreeDistributionModel.EMPIRICAL_TREE_DISTRIBUTION_MODEL,
                new Attribute[] {
                        new Attribute.Default<>(XMLParser.ID, treeModelName),
                        new Attribute.Default<>(EmpiricalTreeDistributionModelParser.FILE_NAME, model.getEmpiricalTreesFilename())
                }, false);

        writer.writeIDref(TaxaParser.TAXA, "taxa"); // @todo - get the actual taxon set for the partition
        writer.writeCloseTag(EmpiricalTreeDistributionModel.EMPIRICAL_TREE_DISTRIBUTION_MODEL);

        writer.writeComment("Statistic to give the current empirical tree");
        writer.writeTag(StatisticParser.STATISTIC,
                new Attribute[]{
                        new Attribute.Default<>(XMLParser.ID, treeModelName + ".currentTree"),
                        new Attribute.Default<>(StatisticParser.NAME,  "Current Tree"),
                }, false);
        writer.writeIDref(EmpiricalTreeDistributionModel.EMPIRICAL_TREE_DISTRIBUTION_MODEL, treeModelName);
        writer.writeCloseTag(StatisticParser.STATISTIC);

        writeTreeModelStatistics(model, writer);
    }


    /**
     * Write default tree model XML block.
     *
     * @param model
     * @param writer the writer
     */
    void writeDefaultTreeModel(PartitionTreeModel model, XMLWriter writer) {

        String prefix = model.getPrefix();

        final String treeModelName = prefix + DefaultTreeModel.TREE_MODEL; // treemodel.treeModel or treeModel

        writer.writeComment("Generate a tree model");
        writer.writeTag(DefaultTreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.ID, treeModelName), false);

        final String STARTING_TREE = InitialTreeGenerator.STARTING_TREE;

        switch (model.getStartingTreeType()) {
            case USER:
                writer.writeIDref("tree", prefix + STARTING_TREE);
                break;
            case UPGMA:
                writer.writeIDref(UPGMATreeParser.UPGMA_TREE, prefix + STARTING_TREE);
                break;
            case RANDOM:
                writer.writeIDref(OldCoalescentSimulatorParser.COALESCENT_TREE, prefix + STARTING_TREE);
                break;
            default:
                throw new IllegalArgumentException("Unknown StartingTreeType");
        }

        writer.writeOpenTag(TreeModelParser.ROOT_HEIGHT);
        writer.writeTag(ParameterParser.PARAMETER,
                new Attribute.Default<String>(XMLParser.ID, treeModelName + "." + TreeModelParser.ROOT_HEIGHT), true);
        writer.writeCloseTag(TreeModelParser.ROOT_HEIGHT);


        writer.writeOpenTag(TreeModelParser.NODE_HEIGHTS, new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"));
        writer.writeTag(ParameterParser.PARAMETER,
                new Attribute.Default<String>(XMLParser.ID, treeModelName + "." + "internalNodeHeights"), true);
        writer.writeCloseTag(TreeModelParser.NODE_HEIGHTS);

        writer.writeOpenTag(TreeModelParser.NODE_HEIGHTS,
                new Attribute[]{
                        new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
                        new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "true")
                });
        writer.writeTag(ParameterParser.PARAMETER,
                new Attribute.Default<String>(XMLParser.ID, treeModelName + "." + "allInternalNodeHeights"), true);
        writer.writeCloseTag(TreeModelParser.NODE_HEIGHTS);

//        int randomLocalClockCount = 0;
//        int autocorrelatedClockCount = 0;
//        for (PartitionData pd : model.getDataPartitions()) { // only the PDs linked to this tree model
//        	PartitionClockModel clockModel = pd.getPartitionClockModel();
//        	switch (clockModel.getClockType()) {
//	        	case AUTOCORRELATED_LOGNORMAL: autocorrelatedClockCount += 1; break;
//	        	case RANDOM_LOCAL_CLOCK: randomLocalClockCount += 1; break;
//        	}
//        }
//
//        if (autocorrelatedClockCount > 1 || randomLocalClockCount > 1 || autocorrelatedClockCount + randomLocalClockCount > 1) {
//        	//FAIL
//            throw new IllegalArgumentException("clock model/tree model combination not implemented by BEAST yet!");
//        }
        // move to validateClockTreeModelCombination(PartitionTreeModel model)

//    	if (autocorrelatedClockCount == 1) {
//        if (count[0] == 1) {
//                writer.writeOpenTag(TreeModelParser.NODE_RATES,
//                        new Attribute[]{
//                                new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "false"),
//                                new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
//                                new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "true")
//                        });
//                writer.writeTag(ParameterParser.PARAMETER,
//                        new Attribute.Default<String>(XMLParser.ID, treeModelName + "." + TreeModelParser.NODE_RATES), true);
//                writer.writeCloseTag(TreeModelParser.NODE_RATES);
//
//                writer.writeOpenTag(TreeModelParser.NODE_RATES,
//                        new Attribute[]{
//                                new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "true"),
//                                new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "false"),
//                                new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "false")
//                        });
//                writer.writeTag(ParameterParser.PARAMETER,
//                        new Attribute.Default<String>(XMLParser.ID,
//                                treeModelName + "." + RateEvolutionLikelihood.ROOTRATE), true);
//                writer.writeCloseTag(TreeModelParser.NODE_RATES);
////    	} else if (randomLocalClockCount == 1 ) {
//        } else

        //+++++++++++++ removed because random local clock XML is changed ++++++++++++++++
//        int[] count = validateClockTreeModelCombination(model);
//        if (count[1] == 1) {
//            writer.writeOpenTag(TreeModelParser.NODE_RATES,
//                    new Attribute[]{
//                            new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "false"),
//                            new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
//                            new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "true")
//                    });
//            writer.writeTag(ParameterParser.PARAMETER,
//                    new Attribute.Default<String>(XMLParser.ID, modelPrefix + ClockType.LOCAL_CLOCK + ".relativeRates"), true);
//            writer.writeCloseTag(TreeModelParser.NODE_RATES);
//
//            writer.writeOpenTag(TreeModelParser.NODE_TRAITS,
//                    new Attribute[]{
//                            new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "false"),
//                            new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
//                            new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "true")
//                    });
//            writer.writeTag(ParameterParser.PARAMETER,
//                    new Attribute.Default<String>(XMLParser.ID, modelPrefix + ClockType.LOCAL_CLOCK + ".changes"), true);
//            writer.writeCloseTag(TreeModelParser.NODE_TRAITS);
//        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_TREE_MODEL, model, writer);

        writer.writeCloseTag(DefaultTreeModel.TREE_MODEL);

        writeTreeModelStatistics(model, writer);
//        if (autocorrelatedClockCount == 1) {
//        if (count[0] == 1) {
//            writer.writeText("");
//            writer.writeOpenTag(CompoundParameter.COMPOUND_PARAMETER,
//                    new Attribute[]{new Attribute.Default<String>(XMLParser.ID, treeModelName + "." + "allRates")});
//            writer.writeIDref(ParameterParser.PARAMETER, treeModelName + "." + TreeModelParser.NODE_RATES);
//            writer.writeIDref(ParameterParser.PARAMETER, treeModelName + "." + RateEvolutionLikelihood.ROOTRATE);
//            writer.writeCloseTag(CompoundParameter.COMPOUND_PARAMETER);
//        }
    }

    void writeTreeModelStatistics(PartitionTreeModel model, XMLWriter writer) {
        String prefix = model.getPrefix();

        final String treeModelName = prefix + DefaultTreeModel.TREE_MODEL; // treemodel.treeModel or treeModel

        writer.writeComment("Statistic for height of the root of the tree");
        writer.writeTag(TreeHeightStatisticParser.TREE_HEIGHT_STATISTIC,
                new Attribute[]{
                        new Attribute.Default<>(XMLParser.ID, prefix + "rootHeight"),
                }, false);
        writer.writeIDref(DefaultTreeModel.TREE_MODEL, treeModelName);
        writer.writeCloseTag(TreeHeightStatisticParser.TREE_HEIGHT_STATISTIC);

        writer.writeComment("Statistic for sum of the branch lengths of the tree (tree length)");
        writer.writeTag(TreeLengthStatisticParser.TREE_LENGTH_STATISTIC,
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, prefix + "treeLength"),
                }, false);
        writer.writeIDref(DefaultTreeModel.TREE_MODEL, treeModelName);
        writer.writeCloseTag(TreeLengthStatisticParser.TREE_LENGTH_STATISTIC);

        if (model.hasTipCalibrations()) {
            writer.writeComment("Statistic for time of most recent common ancestor of tree");
            writer.writeTag(TMRCAStatisticParser.TMRCA_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, prefix + "age(root)"),
                            new Attribute.Default<String>(TMRCAStatisticParser.ABSOLUTE, "true")
                    }, false);
            writer.writeIDref(DefaultTreeModel.TREE_MODEL, treeModelName);
            writer.writeCloseTag(TMRCAStatisticParser.TMRCA_STATISTIC);
        }
    }
}
