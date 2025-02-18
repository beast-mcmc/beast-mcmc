/*
 * TreeLikelihoodGenerator.java
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

import dr.app.beauti.types.ClockType;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodelxml.bigfasttree.thorney.ConstrainedBranchLengthProviderParser;
import dr.evomodelxml.bigfasttree.thorney.PoissonBranchLengthLikelihoodParser;
import dr.evomodelxml.bigfasttree.thorney.ThorneyTreeLikelihoodParser;
import dr.evomodelxml.treedatalikelihood.TreeDataLikelihoodParser;
import dr.evomodelxml.treelikelihood.MarkovJumpsTreeLikelihoodParser;
import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.components.ancestralstates.AncestralStatesComponentOptions;
import dr.app.beauti.options.*;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.bigfasttree.thorney.ConstrainedTreeBranchLengthProvider;
import dr.evomodel.bigfasttree.thorney.ConstrainedTreeModel;
import dr.evomodel.bigfasttree.thorney.MutationBranchMap;
import dr.evomodel.bigfasttree.thorney.PoissonBranchLengthLikelihoodDelegate;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.oldevomodel.sitemodel.GammaSiteModel;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.oldevomodelxml.treelikelihood.AncestralStateTreeLikelihoodParser;
import dr.oldevomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.evoxml.AlignmentParser;
import dr.evoxml.NewickParser;
import dr.evoxml.SimpleTreeParser;
import dr.evoxml.SitePatternsParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.List;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class TreeLikelihoodGenerator extends Generator {

    public TreeLikelihoodGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

    public void writeAllTreeLikelihoods(XMLWriter writer) throws GeneratorException {
        for (List<PartitionData> partitions : options.multiPartitionLists)  {
            writeTreeDataLikelihood(partitions, writer);
            writer.writeText("");
        }

        for (AbstractPartitionData partition : options.otherPartitions) {
            // generate tree likelihoods for the other data partitions
            if (partition.getTaxonList() != null) {
                if (partition instanceof PartitionData) {
                    writeTreeLikelihood((PartitionData) partition, writer);
                    writer.writeText("");

                    // if the partition isn't an instanceof PartitionData then it doesn't
                    // need a TreeLikelihood (it is probably a Tree Partition)
               } else {
                   if(partition instanceof TreePartitionData) {
                        PartitionTreeModel model = partition.getPartitionTreeModel();
                        if(model.isUsingThorneyBEAST()){
                            writeThorneyTreeLikelihood((TreePartitionData) partition, writer);
                            writer.writeText("");
                        }

                   }
                }
            }
        }
    }

    public void writeThorneyTreeLikelihood(TreePartitionData partition,XMLWriter writer){
       // write mutation branch map

/**
 *     <simpleMutationBranchMap id="branchMutationCounts" scale="29903.0">
        <constrainedTreeModel idref="treeModel"/>
        <dataTree>
            <tree idref="dataTree"/>
        </dataTree>
    </simpleMutationBranchMap>
 */
        
        PartitionTreeModel treeModel = partition.getPartitionTreeModel();
        PartitionClockModel clockModel = partition.getPartitionClockModel();
        String prefix = treeModel.getPrefix() + clockModel.getPrefix(); // use the treemodel prefix
        
        String idString = prefix + "treeLikelihood";
        String branchMapIdD = prefix + "branchMutationCounts";


        writer.writeComment("A map between the branches in the tree model and the data providing the number of mutations along each branch.");
        writer.writeTag(ConstrainedBranchLengthProviderParser.MUTATION_BRANCH_MAP, new Attribute[]{
            new Attribute.Default<>(XMLParser.ID,branchMapIdD),
            new Attribute.Default<Double>("scale",treeModel.getThorneyScaler()),
        },false);
            writer.writeIDref(ConstrainedTreeModel.CONSTRAINED_TREE_MODEL,prefix+ DefaultTreeModel.TREE_MODEL);
            writer.writeTag(ConstrainedBranchLengthProviderParser.DATA_TREE, new Attribute[]{}, false);
                writer.writeOpenTag(
                            NewickParser.NEWICK,
                            new Attribute[]{
                                new Attribute.Default<>(XMLParser.ID, prefix + "dataTree"),
                                new Attribute.Default<Boolean>(NewickParser.USING_HEIGHTS, true),
                                new Attribute.Default<Boolean>(NewickParser.USING_DATES, false)
                            }
                    );
                        writer.writeText(TreeUtils.newick(treeModel.getTreePartitionData().getTrees().getTrees().get(0)));
                    writer.writeCloseTag(NewickParser.NEWICK);
            writer.writeCloseTag(ConstrainedBranchLengthProviderParser.DATA_TREE);
        writer.writeCloseTag(ConstrainedBranchLengthProviderParser.MUTATION_BRANCH_MAP);



        writer.writeOpenTag(ThorneyTreeLikelihoodParser.THORNEY_DATA_LIKELIHOOD_DELEGATE,  new Attribute[]{
            new Attribute.Default<String>(XMLParser.ID, idString)
        });

        writer.writeIDref(ConstrainedTreeModel.CONSTRAINED_TREE_MODEL, prefix + DefaultTreeModel.TREE_MODEL);
        writer.writeIDref(ConstrainedBranchLengthProviderParser.MUTATION_BRANCH_MAP, branchMapIdD );
        

        writer.writeTag( PoissonBranchLengthLikelihoodParser.POISSON_BRANCH_LENGTH_LIKELIHOOD,
        new Attribute[] {
            new Attribute.Default<>(XMLParser.ID, prefix + "branchMutationLikelihood"),
            new Attribute.Default<Double>("scale", treeModel.getThorneyScaler()),
    }, true);

        ClockModelGenerator.writeBranchRatesModelRef(clockModel, writer);

        writer.writeCloseTag(ThorneyTreeLikelihoodParser.THORNEY_DATA_LIKELIHOOD_DELEGATE);
    }

    public void writeTreeDataLikelihood(List<PartitionData> partitions, XMLWriter writer) {

        PartitionSubstitutionModel substModel = partitions.get(0).getPartitionSubstitutionModel();
        PartitionTreeModel treeModel = partitions.get(0).getPartitionTreeModel();
        PartitionClockModel clockModel = partitions.get(0).getPartitionClockModel();

        String prefix = treeModel.getPrefix() + clockModel.getPrefix(); // use the treemodel prefix
        String idString = prefix + "treeLikelihood";

        Attribute[] attributes = new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, idString),
                new Attribute.Default<Boolean>(TreeDataLikelihoodParser.USE_AMBIGUITIES, substModel.isUseAmbiguitiesTreeLikelihood()),
                new Attribute.Default<Boolean>(TreeDataLikelihoodParser.USE_PREORDER,
                        clockModel.getClockType() == ClockType.HMC_CLOCK ||
                                clockModel.getClockType() == ClockType.SHRINKAGE_LOCAL_CLOCK)
        };

        writer.writeComment("Likelihood for tree given sequence data");
        writer.writeOpenTag(TreeDataLikelihoodParser.TREE_DATA_LIKELIHOOD, attributes);

        for (PartitionData partition : partitions) {
            substModel = partition.getPartitionSubstitutionModel();
            PartitionClockModel cm = partition.getPartitionClockModel();
            if (clockModel == null) {
                clockModel = cm;
            } else if (clockModel != cm) {
                throw new RuntimeException("All the partitions in a TreeDataLikelihood should share the same clock model.");
            }

            if (substModel.getCodonHeteroPattern() != null) {

                for (int num = 1; num <= substModel.getCodonPartitionCount(); num++) {
                    String prefix1 = partition.getPrefix() + substModel.getPrefixCodon(num);

                    writer.writeOpenTag(TreeDataLikelihoodParser.PARTITION);
                    writeCodonPatternsRef(prefix1, num, substModel.getCodonPartitionCount(), writer);

                    writer.writeIDref(GammaSiteModel.SITE_MODEL, substModel.getPrefix(num) + SiteModel.SITE_MODEL);

                    writer.writeCloseTag(TreeDataLikelihoodParser.PARTITION);
                }


            } else {
                String prefix1 = partition.getPrefix();
                writer.writeOpenTag(TreeDataLikelihoodParser.PARTITION);
                writer.writeIDref(SitePatternsParser.PATTERNS, prefix1 + SitePatternsParser.PATTERNS);
                writer.writeIDref(GammaSiteModel.SITE_MODEL, substModel.getPrefix() + SiteModel.SITE_MODEL);
                writer.writeCloseTag(TreeDataLikelihoodParser.PARTITION);
            }


        }

        writer.writeIDref(DefaultTreeModel.TREE_MODEL, treeModel.getPrefix() + DefaultTreeModel.TREE_MODEL);
        ClockModelGenerator.writeBranchRatesModelRef(clockModel, writer);


        writer.writeCloseTag(TreeDataLikelihoodParser.TREE_DATA_LIKELIHOOD);

    }

    /**
     * Write the tree likelihood XML block.
     *
     * @param partition the partition  to write likelihood block for
     * @param writer    the writer
     */
    public void writeTreeLikelihood(PartitionData partition, XMLWriter writer) {

        AncestralStatesComponentOptions ancestralStatesOptions = (AncestralStatesComponentOptions) options
                .getComponentOptions(AncestralStatesComponentOptions.class);

        PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();

        if (model.isDolloModel()) {
            return; // DolloComponent will add tree likelihood
        }

        String treeLikelihoodTag = TreeLikelihoodParser.TREE_LIKELIHOOD;
        if (ancestralStatesOptions.usingAncestralStates(partition)) {
            treeLikelihoodTag = TreeLikelihoodParser.ANCESTRAL_TREE_LIKELIHOOD;
            if (ancestralStatesOptions.isCountingStates(partition)) {
                // State change counting uses the MarkovJumpsTreeLikelihood but
                // dNdS robust counting doesn't as it has its own counting code...
                if (!ancestralStatesOptions.dNdSRobustCounting(partition)) {
                    treeLikelihoodTag = MarkovJumpsTreeLikelihoodParser.MARKOV_JUMP_TREE_LIKELIHOOD;
                }
            }
        }

        if (model.getDataType().getType() == DataType.NUCLEOTIDES && model.getCodonHeteroPattern() != null) {

            for (int i = 1; i <= model.getCodonPartitionCount(); i++) {
                writeTreeLikelihood(treeLikelihoodTag, TreeLikelihoodParser.TREE_LIKELIHOOD, i, partition, writer);
            }


        } else {
            writeTreeLikelihood(treeLikelihoodTag, TreeLikelihoodParser.TREE_LIKELIHOOD, -1, partition, writer);
        }
    }

    /**
     * Write the tree likelihood XML block.
     *
     * @param id        the id of the tree likelihood
     * @param num       the likelihood number
     * @param partition the partition to write likelihood block for
     * @param writer    the writer
     */
    private void writeTreeLikelihood(String tag, String id, int num, PartitionData partition, XMLWriter writer) {

        PartitionSubstitutionModel substModel = partition.getPartitionSubstitutionModel();
        PartitionTreeModel treeModel = partition.getPartitionTreeModel();
        PartitionClockModel clockModel = partition.getPartitionClockModel();

        writer.writeComment("Likelihood for tree given sequence data");

        String prefix;
        if (num > 0) {
            prefix = partition.getPrefix() + substModel.getPrefixCodon(num);
        } else {
            prefix = partition.getPrefix();
        }

        String idString = prefix + id;

        Attribute[] attributes;
        if (tag.equals(MarkovJumpsTreeLikelihoodParser.MARKOV_JUMP_TREE_LIKELIHOOD)) {
            AncestralStatesComponentOptions ancestralStatesOptions = (AncestralStatesComponentOptions) options
                    .getComponentOptions(AncestralStatesComponentOptions.class);
            boolean saveCompleteHistory = ancestralStatesOptions.isCompleteHistoryLogging(partition);
            attributes = new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, idString),
                    new Attribute.Default<Boolean>(TreeLikelihoodParser.USE_AMBIGUITIES, substModel.isUseAmbiguitiesTreeLikelihood()),
                    new Attribute.Default<Boolean>(MarkovJumpsTreeLikelihoodParser.USE_UNIFORMIZATION, true),
                    new Attribute.Default<Integer>(MarkovJumpsTreeLikelihoodParser.NUMBER_OF_SIMULANTS, 1),
                    new Attribute.Default<String>(AncestralStateTreeLikelihoodParser.RECONSTRUCTION_TAG_NAME, prefix + AncestralStateTreeLikelihoodParser.RECONSTRUCTION_TAG),
                    new Attribute.Default<String>(MarkovJumpsTreeLikelihoodParser.SAVE_HISTORY, saveCompleteHistory ? "true" : "false"),
                    new Attribute.Default<String>(MarkovJumpsTreeLikelihoodParser.LOG_HISTORY, saveCompleteHistory ? "true" : "false"),
            };
        } else if (tag.equals(TreeLikelihoodParser.ANCESTRAL_TREE_LIKELIHOOD)) {
            attributes = new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, idString),
                    new Attribute.Default<Boolean>(TreeLikelihoodParser.USE_AMBIGUITIES, substModel.isUseAmbiguitiesTreeLikelihood()),
                    new Attribute.Default<String>(AncestralStateTreeLikelihoodParser.RECONSTRUCTION_TAG_NAME, prefix + AncestralStateTreeLikelihoodParser.RECONSTRUCTION_TAG),
            };
        } else {
            attributes = new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, idString),
                    new Attribute.Default<Boolean>(TreeLikelihoodParser.USE_AMBIGUITIES, substModel.isUseAmbiguitiesTreeLikelihood())
            };
        }

        writer.writeOpenTag(tag, attributes);

        if (!options.samplePriorOnly) {
            if (num > 0) {
                writeCodonPatternsRef(prefix, num, substModel.getCodonPartitionCount(), writer);
            } else {
                writer.writeIDref(SitePatternsParser.PATTERNS, prefix + SitePatternsParser.PATTERNS);
            }
        } else {
            // We just need to use the dummy alignment
            writer.writeIDref(AlignmentParser.ALIGNMENT, partition.getAlignment().getId());
        }

        writer.writeIDref(DefaultTreeModel.TREE_MODEL, treeModel.getPrefix() + DefaultTreeModel.TREE_MODEL);

        if (num > 0) {
            writer.writeIDref(GammaSiteModel.SITE_MODEL, substModel.getPrefix(num) + SiteModel.SITE_MODEL);
        } else {
            writer.writeIDref(GammaSiteModel.SITE_MODEL, substModel.getPrefix() + SiteModel.SITE_MODEL);
        }

        ClockModelGenerator.writeBranchRatesModelRef(clockModel, writer);

        generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_TREE_LIKELIHOOD, partition, prefix, writer);

        writer.writeCloseTag(tag);
    }

    /**
     * Write a list of idrefs to all treelikelihoods
     * @param writer
     */
    public void writeTreeLikelihoodReferences(XMLWriter writer) {
        for (List<PartitionData> partitions : options.multiPartitionLists)  {
            // TreeDataLikelihoods are labelled by their tree and clock models
            PartitionTreeModel treeModel = partitions.get(0).getPartitionTreeModel();
            PartitionClockModel clockModel = partitions.get(0).getPartitionClockModel();

            String prefix = treeModel.getPrefix() + clockModel.getPrefix(); // use the treemodel prefix
            String idString = prefix + "treeLikelihood";

            writer.writeIDref(TreeDataLikelihoodParser.TREE_DATA_LIKELIHOOD, idString);
        }

        for (AbstractPartitionData partition : options.otherPartitions) {
            // generate tree likelihoods for the other data partitions

            // TreeLikelihood are labelled by their partition
            AncestralStatesComponentOptions ancestralStatesOptions = (AncestralStatesComponentOptions) options
                    .getComponentOptions(AncestralStatesComponentOptions.class);

            // TODO moving poisson branch length to a substitution model would make this logic simplier

            PartitionTreeModel treeModel = partition.getPartitionTreeModel();
            if(treeModel.isUsingThorneyBEAST()){
                PartitionClockModel clockModel = partition.getPartitionClockModel();
                String prefix = treeModel.getPrefix() + clockModel.getPrefix(); // use the treemodel prefix
                String idString = prefix + "treeLikelihood"; // Line 117
                writer.writeIDref(ThorneyTreeLikelihoodParser.THORNEY_DATA_LIKELIHOOD_DELEGATE, idString);
            
            }
            PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();

            if (model == null) {
                return; // some partitions don't have evolutionary models (i.e., empirical trees)
            }

            if (model.isDolloModel()) {
                return; // DolloComponent will add tree likelihood
            }

            String treeLikelihoodTag = TreeLikelihoodParser.TREE_LIKELIHOOD;
            if (ancestralStatesOptions.usingAncestralStates(partition)) {
                treeLikelihoodTag = TreeLikelihoodParser.ANCESTRAL_TREE_LIKELIHOOD;
                if (ancestralStatesOptions.isCountingStates(partition)) {
                    // State change counting uses the MarkovJumpsTreeLikelihood but
                    // dNdS robust counting doesn't as it has its own counting code...
                    if (!ancestralStatesOptions.dNdSRobustCounting(partition)) {
                        treeLikelihoodTag = MarkovJumpsTreeLikelihoodParser.MARKOV_JUMP_TREE_LIKELIHOOD;
                    }
                }
            }

            if (model.getDataType().getType() == DataType.NUCLEOTIDES && model.getCodonHeteroPattern() != null) {

                for (int i = 1; i <= model.getCodonPartitionCount(); i++) {
                    writeTreeLikelihoodRef(treeLikelihoodTag, TreeLikelihoodParser.TREE_LIKELIHOOD, i, partition, writer);
                }


            } else {
                writeTreeLikelihoodRef(treeLikelihoodTag, TreeLikelihoodParser.TREE_LIKELIHOOD, -1, partition, writer);
            }
        }
    }

    private void writeTreeLikelihoodRef(String tag, String id, int num, AbstractPartitionData partition, XMLWriter writer) {

        PartitionSubstitutionModel substModel = partition.getPartitionSubstitutionModel();

        String prefix;
        if (num > 0) {
            prefix = partition.getPrefix() + substModel.getPrefixCodon(num);
        } else {
            prefix = partition.getPrefix();
        }

        String idString = prefix + id;

        writer.writeIDref(tag, idString);
    }

    /**
     * Is the data partition suitable for using the MultipartitionTreeDataLikelihoodDelegate.
     *
     * @param partition the partition  to write likelihood block for
     */
    public boolean canUseMultiPartition(AbstractPartitionData partition) {

        AncestralStatesComponentOptions ancestralStatesOptions = (AncestralStatesComponentOptions) options
                .getComponentOptions(AncestralStatesComponentOptions.class);

        PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();

        if (model == null || model.isDolloModel()) {
            return false;
        }

        if (ancestralStatesOptions.usingAncestralStates(partition)) {
            return false;
        }

        if (model.getDataType().getType() != DataType.NUCLEOTIDES && model.getDataType().getType() != DataType.AMINO_ACIDS) {
            return false;
        }

        return true;
    }

}
