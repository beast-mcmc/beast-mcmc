/*
 * PartitionTreeModel.java
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

package dr.app.beauti.options;

import dr.app.beauti.types.*;
import dr.evolution.tree.Tree;

import java.util.Arrays;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class PartitionTreeModel extends PartitionOptions {

    private static final long serialVersionUID = 4829401415152235625L;
    private static final String DEFAULT_EMPIRICAL_TREES_FILENAME = "empirical.trees";

    private PartitionTreePrior treePrior;

    private StartingTreeType startingTreeType = StartingTreeType.RANDOM;
    private Tree userStartingTree = null;

    private TreeAsDataType treeAsDataType = TreeAsDataType.EMPRICAL_TREES;

    private boolean isUsingExternalEmpiricalTreeFile = false;
    private String empiricalTreesFilename = null;

    private double thorneyScaler = 1.0;
    private boolean isNewick = true;

    private boolean hasTipCalibrations = false;
    private boolean hasNodeCalibrations = false;
    public static List<ClockType> thorneyCompatibleClocks = Arrays.asList( ClockType.STRICT_CLOCK, ClockType.UNCORRELATED,
            ClockType.FIXED_LOCAL_CLOCK, ClockType.RANDOM_LOCAL_CLOCK, ClockType.FIXED_LOCAL_CLOCK,
            ClockType.MIXED_EFFECTS_CLOCK );

    private final TreePartitionData treePartitionData;

    public PartitionTreeModel(BeautiOptions options, String name) {
        super(options, name);

        treePartitionData = null;

        initModelParametersAndOpererators();
    }

    public TreePartitionData getTreePartitionData() {
        return treePartitionData;
    }

    public PartitionTreeModel(BeautiOptions options, TreePartitionData treePartitionData) {
        super(options, treePartitionData.getName());

        this.treePartitionData = treePartitionData;
        initModelParametersAndOpererators();
    }

    /**
     * A copy constructor
     *
     * @param options the beauti options
     * @param name    the name of the new model
     * @param source  the source model
     */
    public PartitionTreeModel(BeautiOptions options, String name, PartitionTreeModel source) {
        super(options, name);

        treePartitionData = null;

        treePrior = source.treePrior;
        startingTreeType = source.startingTreeType;
        userStartingTree = source.userStartingTree;

        isNewick = source.isNewick;
        // initialRootHeight = source.initialRootHeight;

        initModelParametersAndOpererators();
    }

    public void initModelParametersAndOpererators() {

        createParameter("tree", "The tree");
        createParameter("treeModel.internalNodeHeights", "internal node heights of the tree (except the root)");
        createParameter("treeModel.allInternalNodeHeights", "internal node heights of the tree");
        createParameterTree(this, "treeModel.rootHeight", "root height of the tree", true);

        if (treePartitionData == null) {
            // TODO treeBitMove should move to PartitionClockModelTreeModelLink, after
            // Alexei finish
            createOperator("treeBitMove", "Tree", "Swaps the rates and change locations of local clocks", "tree",
                    OperatorType.TREE_BIT_MOVE, -1.0, treeWeights);

            createScaleOperator("treeModel.rootHeight", demoTuning, demoWeights);
            createOperator("uniformHeights", "Internal node heights", "Draws new internal node heights uniformally",
                    "treeModel.internalNodeHeights", OperatorType.UNIFORM, -1, branchWeights);

            // This scale operator is used instead of the up/down if the rate is fixed.
            new Operator.Builder("treeModel.allInternalNodeHeights", "Scales all internal node heights in tree",
                    getParameter("treeModel.allInternalNodeHeights"), OperatorType.SCALE_ALL, 0.75, rateWeights)
                    .build(operators);

            createOperator("subtreeSlide", "Tree", "Performs the subtree-slide rearrangement of the tree", "tree",
                    OperatorType.SUBTREE_SLIDE, 1.0, treeWeights);
            createOperator("narrowExchange", "Tree", "Performs local rearrangements of the tree", "tree",
                    OperatorType.NARROW_EXCHANGE, -1, treeWeights);
            createOperator("wideExchange", "Tree", "Performs global rearrangements of the tree", "tree",
                    OperatorType.WIDE_EXCHANGE, -1, demoWeights);
            createOperator("wilsonBalding", "Tree", "Performs the Wilson-Balding rearrangement of the tree", "tree",
                    OperatorType.WILSON_BALDING, -1, demoWeights);

            double weight = Math.max(options.taxonList.getTaxonCount(), 30);
            createOperator("subtreeLeap", "Tree", "Performs the subtree-leap rearrangement of the tree", "tree",
                    OperatorType.SUBTREE_LEAP, 1.0, weight);

            weight = Math.max(weight / 10, 3);
            createOperator("FHSPR", "Tree", "Performs the fixed-height subtree prune/regraft of the tree", "tree",
                    OperatorType.FIXED_HEIGHT_SUBTREE_PRUNE_REGRAFT, 1.0, weight);
        } else {
            // Thorney BEAST operators
            // Big weights operators
            double weight = Math.max(options.taxonList.getTaxonCount(), 30);

            // does STL works on constrained trees?
            createOperator("subtreeLeap", "Tree", "Performs the subtree-leap rearrangement of the tree", "tree",
                    OperatorType.SUBTREE_LEAP, 1.0, weight);
            createOperator("uniformSPG", "Tree",
                    "Performs the subtree prune regraft within a polytomony rearrangement of the tree", "tree",
                    OperatorType.UNIFORM_SUBTREE_PRUNE_REGRAFT, -1, weight);

            createOperator("treeBitMove", "Tree", "Swaps the rates and change locations of local clocks", "tree",
                    OperatorType.TREE_BIT_MOVE, -1.0, treeWeights);

            createOperator("uniformHeights", "Internal node heights", "Draws new internal node heights uniformally",
                    "treeModel.internalNodeHeights", OperatorType.UNIFORM, -1, weight);

            // This scale operator is used instead of the up/down if the rate is fixed.
            new Operator.Builder("treeModel.allInternalNodeHeights", "Scales all internal node heights in tree",
                    getParameter("treeModel.allInternalNodeHeights"), OperatorType.SCALE_ALL, 0.75, rateWeights)
                    .build(operators);

            // In the TB setting this will almost always give the same likelihood since
            // coalescent doesn't change
            // and if we are in the non mutation zone the treelikelihood will be the same.
            // not compatible with constrained tree right now.
            // createOperator("FHSPR", "Tree", "Performs the fixed-height subtree
            // prune/regraft of the tree", "tree",
            // OperatorType.FIXED_HEIGHT_SUBTREE_PRUNE_REGRAFT, 1.0, weight);

            createOperator("narrowExchange", "Tree", "Performs local rearrangements of the tree", "tree",
                    OperatorType.NARROW_EXCHANGE, -1, weight);

            weight = Math.max(weight / 100, 3);
            createScaleOperator("treeModel.rootHeight", demoTuning, weight);

            createOperator("wideExchange", "Tree", "Performs global rearrangements of the tree", "tree",
                    OperatorType.WIDE_EXCHANGE, -1, weight);

            createOperator("wilsonBalding", "Tree", "Performs the Wilson-Balding rearrangement of the tree", "tree",
                    OperatorType.WILSON_BALDING, -1, weight);
        }

        createOperator("empiricalTreeSwap", "Tree", "Sets the current tree from the empirical set", "tree",
                OperatorType.EMPIRICAL_TREE_SWAP, -1.0, treeWeights);
    }

    @Override
    public List<Parameter> selectParameters(List<Parameter> parameters) {
        // setAvgRootAndRate();

        // Don't add these to the parameter list (as they don't appear in the table),
        // but call
        // get parameter so their id prefix can be set.

        getParameter("tree");
        getParameter("treeModel.internalNodeHeights");
        getParameter("treeModel.allInternalNodeHeights");

        Parameter rootHeightParameter = getParameter("treeModel.rootHeight");
        if (rootHeightParameter.priorType == PriorType.NONE_TREE_PRIOR || !rootHeightParameter.isPriorEdited()) {
            rootHeightParameter.setInitial(getInitialRootHeight());
            rootHeightParameter.truncationLower = options.maximumTipHeight;
            rootHeightParameter.uniformLower = options.maximumTipHeight;
            rootHeightParameter.isTruncated = true;
        }

        rootHeightParameter.isCalibratedYule = treePrior.getNodeHeightPrior() == TreePriorType.YULE_CALIBRATION;
        parameters.add(rootHeightParameter);

        return parameters;
    }

    @Override
    public List<Operator> selectOperators(List<Operator> operators) {
        // setAvgRootAndRate();

        if (isUsingEmpiricalTrees()) {
            operators.add(getOperator("empiricalTreeSwap"));
        } else if (isUsingThorneyBEAST()) {

            operators.add(getOperator("treeModel.rootHeight"));
            operators.add(getOperator("uniformHeights"));
            // TODO check if there are polytomies otherwise these will error out
            operators.add(getOperator("uniformSPG"));
            operators.add(getOperator("narrowExchange"));
            operators.add(getOperator("wideExchange"));
            operators.add(getOperator("wilsonBalding"));
            // TODO add back in once operators is functional
            // operators.add(getOperator("FHSPR"));

        } else {
            if (treePartitionData == null) {
                Operator subtreeSlideOp = getOperator("subtreeSlide");
                if (!subtreeSlideOp.isTuningEdited()) {
                    double tuning = 1.0;
                    if (!Double.isNaN(getInitialRootHeight()) && !Double.isInfinite(getInitialRootHeight())) {
                        tuning = getInitialRootHeight() / 10.0;
                    }
                    subtreeSlideOp.setTuning(tuning);
                }

                operators.add(subtreeSlideOp);
                operators.add(getOperator("narrowExchange"));
                operators.add(getOperator("wideExchange"));
                operators.add(getOperator("wilsonBalding"));

                operators.add(getOperator("treeModel.rootHeight"));
                operators.add(getOperator("uniformHeights"));

                operators.add(getOperator("subtreeLeap"));
                operators.add(getOperator("FHSPR"));

                if (options.operatorSetType != OperatorSetType.CUSTOM) {
                    // do nothing
                    boolean defaultInUse = false;
                    boolean branchesInUse = false;
                    boolean newTreeOperatorsInUse = false;
                    boolean adaptiveMultivariateInUse = false;

                    // if not a fixed tree then sample tree space
                    if (options.operatorSetType != OperatorSetType.FIXED_TREE) {
                        if (options.operatorSetType == OperatorSetType.DEFAULT) {
                            newTreeOperatorsInUse = true; // default is now the new tree operators
                        } else if (options.operatorSetType == OperatorSetType.CLASSIC) {
                            defaultInUse = true;
                            branchesInUse = true;
                        } else if (options.operatorSetType == OperatorSetType.FIXED_TREE_TOPOLOGY) {
                            branchesInUse = true;
                        } else if (options.operatorSetType == OperatorSetType.ADAPTIVE_MULTIVARIATE) {
                            newTreeOperatorsInUse = true;
                            adaptiveMultivariateInUse = true;
                        } else {
                            throw new IllegalArgumentException("Unknown operator set type");
                        }
                    }

                    getOperator("subtreeSlide").setUsed(defaultInUse);
                    getOperator("narrowExchange").setUsed(defaultInUse);
                    getOperator("wideExchange").setUsed(defaultInUse);
                    getOperator("wilsonBalding").setUsed(defaultInUse);

                    getOperator("treeModel.rootHeight").setUsed(branchesInUse);
                    getOperator("treeModel.allInternalNodeHeights").setUsed(branchesInUse);
                    getOperator("uniformHeights").setUsed(branchesInUse);

                    getOperator("subtreeLeap").setUsed(newTreeOperatorsInUse);
                    getOperator("FHSPR").setUsed(newTreeOperatorsInUse);
                }
            } else {
                getOperator("subtreeLeap").setUsed(true);
                getOperator("FHSPR").setUsed(true);
            }
        }
        return operators;
    }

    /////////////////////////////////////////////////////////////

    public PartitionTreePrior getPartitionTreePrior() {
        return treePrior;
    }

    public void setPartitionTreePrior(PartitionTreePrior treePrior) {
        options.clearDataPartitionCaches();
        this.treePrior = treePrior;
    }

    public StartingTreeType getStartingTreeType() {
        return startingTreeType;
    }

    public void setStartingTreeType(StartingTreeType startingTreeType) {
        this.startingTreeType = startingTreeType;
    }

    public Tree getUserStartingTree() {
        return userStartingTree;
    }

    public void setUserStartingTree(Tree userStartingTree) {
        this.userStartingTree = userStartingTree;
    }

    public void setTreeAsDataType(TreeAsDataType treeAsDataType) {
        this.treeAsDataType = treeAsDataType;
    }

    public TreeAsDataType getTreeAsDataType() {
        return treeAsDataType;
    }

    public boolean isUsingEmpiricalTrees() {
        return treePartitionData != null && treeAsDataType == TreeAsDataType.EMPRICAL_TREES;
    }

    public boolean isUsingThorneyBEAST() {
        return treePartitionData != null && treeAsDataType == TreeAsDataType.THORNEY_BEAST;
    }

    public void setUsingExternalEmpiricalTreeFile(boolean isUsingExternalEmpiricalTreeFile) {
        this.isUsingExternalEmpiricalTreeFile = isUsingExternalEmpiricalTreeFile;
    }

    public boolean isUsingExternalEmpiricalTreeFile() {
        return isUsingExternalEmpiricalTreeFile;
    }

    public String getEmpiricalTreesFilename() {
        if (empiricalTreesFilename == null || empiricalTreesFilename.isEmpty()) {
            empiricalTreesFilename = DEFAULT_EMPIRICAL_TREES_FILENAME;
        }
        return empiricalTreesFilename;
    }

    public void setEmpiricalTreesFilename(String empiricalTreesFilename) {
        this.empiricalTreesFilename = empiricalTreesFilename;
    }

    public double getThorneyScaler() {
        return thorneyScaler;
    }

    public void setThorneyScaler(double scaler) {
        this.thorneyScaler = scaler;
    }

    public void setTipCalibrations(boolean hasTipCalibrations) {
        this.hasTipCalibrations = hasTipCalibrations;
    }

    public boolean hasTipCalibrations() {
        return hasTipCalibrations;
    }

    public void setNodeCalibrations(boolean hasNodeCalibrations) {
        this.hasNodeCalibrations = hasNodeCalibrations;
    }

    public boolean hasNodeCalibrations() {
        return hasNodeCalibrations;
    }

    public double getInitialRootHeight() {
        return Double.NaN;
        // return getAvgRootAndRate()[0];
    }

    // public void setInitialRootHeight(double initialRootHeight) {
    // this.initialRootHeight = initialRootHeight;
    // }

    // private void calculateInitialRootHeightPerTree() {
    // initialRootHeight = options.clockModelOptions
    // .calculateInitialRootHeightAndRate(options.getDataPartitions(this)) [0];
    // }

    public String getPrefix() {
        String prefix = "";
        if (options.getPartitionTreeModels().size() > 1) { // || options.isSpeciesAnalysis()
            // There is more than one active partition model
            prefix += getName() + ".";
        }
        return prefix;
    }

    public int getDimension() { // n-1
        return options.getTaxonCount(options.getDataPartitions(this)) - 1;
    }

    public int getTaxonCount() {
        return options.getTaxonCount(options.getDataPartitions(this));
    }
}
