/*
 * OperatorsGenerator.java
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
import dr.app.beauti.options.*;
import dr.app.beauti.types.ClockType;
import dr.app.beauti.types.TreePriorType;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.datatype.DataType;
import dr.evomodel.branchratemodel.BranchSpecificFixedEffects;
import dr.evomodel.bigfasttree.thorney.ConstrainedTreeModel;
import dr.evomodel.operators.BitFlipInSubstitutionModelOperator;
import dr.evomodel.operators.EmpiricalTreeDistributionOperator;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.EmpiricalTreeDistributionModel;
import dr.evomodelxml.bigfasttree.thorney.UniformSubtreePruneRegraftParser;
import dr.evomodelxml.branchratemodel.AutoCorrelatedBranchRatesDistributionParser;
import dr.evomodelxml.branchratemodel.AutoCorrelatedGradientWrtIncrementsParser;
import dr.evomodelxml.branchratemodel.BranchRateGradientWrtIncrementsParser;
import dr.evomodelxml.coalescent.GMRFSkyrideLikelihoodParser;
import dr.evomodelxml.coalescent.operators.GMRFSkyrideBlockUpdateOperatorParser;
import dr.evomodelxml.coalescent.operators.SampleNonActiveGibbsOperatorParser;
import dr.evomodelxml.continuous.hmc.BranchRateGradientParser;
import dr.evomodelxml.continuous.hmc.LocationScaleGradientParser;
import dr.evomodelxml.operators.*;
import dr.evomodelxml.treedatalikelihood.TreeDataLikelihoodParser;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.ParameterParser;
import dr.inference.operators.*;
import dr.inferencexml.SignTransformParser;
import dr.inferencexml.distribution.shrinkage.BayesianBridgeDistributionModelParser;
import dr.inferencexml.hmc.CompoundGradientParser;
import dr.inferencexml.hmc.HessianWrapperParser;
import dr.inferencexml.hmc.JointGradientParser;
import dr.inferencexml.model.CompoundParameterParser;
import dr.inferencexml.operators.*;
import dr.inferencexml.operators.hmc.HamiltonianMonteCarloOperatorParser;
import dr.inferencexml.operators.shrinkage.BayesianBridgeShrinkageOperatorParser;
import dr.oldevomodel.substmodel.AbstractSubstitutionModel;
import dr.oldevomodelxml.substmodel.GeneralSubstitutionModelParser;
import dr.util.Attribute;
import dr.util.Transform;
import dr.util.TransformParsers;
import dr.xml.XMLParser;

import java.util.ArrayList;
import java.util.List;

import static dr.inferencexml.distribution.PriorParsers.GAMMA_PRIOR;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class OperatorsGenerator extends Generator {

    public static final boolean NEW_AVMVN = true;

    public OperatorsGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

    /**
     * Write the operator schedule XML block.
     *
     * @param operators the list of operators
     * @param writer    the writer
     */
    public void writeOperatorSchedule(List<Operator> operators, XMLWriter writer) {
        Attribute[] operatorAttributes;

        // certain models would benefit from a logarithm operator optimization
        boolean shouldLogCool = false;
        for (PartitionTreePrior partition : options.getPartitionTreePriors()) {
            if (partition.getNodeHeightPrior() == TreePriorType.SKYGRID ||
                    partition.getNodeHeightPrior() == TreePriorType.GMRF_SKYRIDE) {
                shouldLogCool = true;
                break;
            }
        }
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
            if (model.getDataType().getType() == DataType.GENERAL ||
                    model.getDataType().getType() == DataType.CONTINUOUS) {
                shouldLogCool = true;
                break;
            }
        }

        operatorAttributes = new Attribute[] {
                new Attribute.Default<String>(XMLParser.ID, "operators"),
                new Attribute.Default<String>(SimpleOperatorScheduleParser.OPTIMIZATION_SCHEDULE,
                        (shouldLogCool ?
                                OperatorSchedule.OptimizationTransform.LOG.toString() :
                                OperatorSchedule.DEFAULT_TRANSFORM.toString()))
        };

        writer.writeComment("Define operators");
        writer.writeOpenTag(
                SimpleOperatorScheduleParser.OPERATOR_SCHEDULE,
                operatorAttributes
//				new Attribute[]{new Attribute.Default<String>(XMLParser.ID, "operators")}
        );

        for (Operator operator : operators) {
            if (operator.getWeight() > 0. && operator.isUsed()) {
                writeOperator(operator, writer);
            }
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_OPERATORS, writer); // Added for special operators

        writer.writeCloseTag(SimpleOperatorScheduleParser.OPERATOR_SCHEDULE);
    }

    private void writeOperator(Operator operator, XMLWriter writer) {

        String prefix = operator.getPrefix();

        switch (operator.getOperatorType()) {

            case SCALE:
                writeScaleOperator(operator, writer);
                break;
            case RANDOM_WALK:
                writeRandomWalkOperator(operator, writer);
                break;
            case RANDOM_WALK_ABSORBING:
                writeRandomWalkOperator(operator, RandomWalkOperator.BoundaryCondition.absorbing, writer);
                break;
            case RANDOM_WALK_REFLECTING:
                writeRandomWalkOperator(operator, RandomWalkOperator.BoundaryCondition.reflecting, writer);
                break;
            case RANDOM_WALK_LOG:
                writeRandomWalkOperator(operator, RandomWalkOperator.BoundaryCondition.log, writer);
                break;
            case RANDOM_WALK_LOGIT:
                writeRandomWalkOperator(operator, RandomWalkOperator.BoundaryCondition.logit, writer);
                break;
            case INTEGER_RANDOM_WALK:
                writeIntegerRandomWalkOperator(operator, writer);
                break;
            case UP_DOWN:
                writeUpDownOperator(UpDownOperatorParser.UP_DOWN_OPERATOR, operator, writer);
                break;
            case SCALE_ALL:
                writeScaleAllOperator(operator, writer);
                break;
            case SCALE_INDEPENDENTLY:
                writeScaleOperator(operator, writer, true);
                break;
            case DELTA_EXCHANGE:
                writeDeltaOperator(operator, false, writer);
                break;
            case WEIGHTED_DELTA_EXCHANGE:
                writeDeltaOperator(operator, true, writer);
                break;
            case INTEGER_DELTA_EXCHANGE:
                writeIntegerDeltaOperator(operator, writer);
                break;
            case SWAP:
                writeSwapOperator(operator, writer);
                break;
            case BITFLIP:
                writeBitFlipOperator(operator, writer);
                break;
            case BITFLIP_IN_SUBST:
                writeBitFlipInSubstOperator(operator, writer);
                break;
            case RATE_BIT_EXCHANGE:
                writeRateBitExchangeOperator(operator, writer);
                break;
            case TREE_BIT_MOVE:
                writeTreeBitMoveOperator(operator, prefix, writer);
                break;
            case UNIFORM:
                writeUniformOperator(operator, writer);
                break;
            case INTEGER_UNIFORM:
                writeIntegerUniformOperator(operator, writer);
                break;
            case SUBTREE_LEAP:
                writeSubtreeLeapOperator(operator, prefix, writer);
                break;
            case FIXED_HEIGHT_SUBTREE_PRUNE_REGRAFT:
                writeFHSPROperator(operator, prefix, writer);
                break;
            case EMPIRICAL_TREE_SWAP:
                writeEmpiricalTreeSwapOperator(operator, prefix, writer);
                break;
            case SUBTREE_SLIDE:
                writeSubtreeSlideOperator(operator, prefix, writer);
                break;
            // write multivariate operator
            case NARROW_EXCHANGE:
                writeNarrowExchangeOperator(operator, prefix, writer);
                break;
            case WIDE_EXCHANGE:
                writeWideExchangeOperator(operator, prefix, writer);
                break;
            case WILSON_BALDING:
                writeWilsonBaldingOperator(operator, prefix, writer);
                break;
            case SAMPLE_NONACTIVE:
                writeSampleNonActiveOperator(operator, writer);
                break;
            case SCALE_WITH_INDICATORS:
                writeScaleWithIndicatorsOperator(operator, writer);
                break;
            case GMRF_BLOCKUPDATE_OPERATOR:
                writeGMRFBlockUpdateOperator(operator, prefix, writer);
                break;
            case SKY_GRID_BLOCKUPDATE_OPERATOR:
                writeSkyGridBlockUpdateOperator(operator, prefix, writer);
                break;
            case SKY_GRID_HMC_OPERATOR:
                writeSkyGridHMCOperator(operator, prefix, writer);
                break;
            case ADAPTIVE_MULTIVARIATE:
                writeAdaptiveMultivariateOperator(operator, writer);
                break;
            case RELAXED_CLOCK_HMC_RATE_OPERATOR:
                writeRelaxedClockHMCRateOperator(operator, prefix,writer);
                break;
            case RELAXED_CLOCK_HMC_SCALE_OPERATOR:
                writeRelaxedClockHMCScaleOperator(operator, prefix,writer);
                break;
            case SHRINKAGE_CLOCK_HMC_OPERATOR:
                writeShrinkageClockHMCOperator(operator, prefix, writer);
                break;
            case SHRINKAGE_CLOCK_GIBBS_OPERATOR:
                writeShrinkageClockGibbsOperator(operator, prefix, writer);
                break;
            case NODE_HEIGHT_OPERATOR_UNIFORM :
                writeBFTUniformNodeHeightOperator(operator, prefix, writer);
                break;
            case NODE_HEIGHT_OPERATOR_ROOT:
                writeBFTRootScaleOperator(operator, prefix, writer);
                break;
            case UNIFORM_SUBTREE_PRUNE_REGRAFT:
                writeUniformSPROperator(operator, prefix, writer);
                break;
            default:
                throw new IllegalArgumentException("Unknown operator type");
        }
    }

    private void writeParameter1Ref(XMLWriter writer, Operator operator) {
        writer.writeIDref(ParameterParser.PARAMETER, operator.getParameter1().getName());
    }

    private void writeParameter2Ref(XMLWriter writer, Operator operator) {
        writer.writeIDref(ParameterParser.PARAMETER, operator.getParameter2().getName());
    }

    private void writeOperatorRef(XMLWriter writer, Operator operator) {
        writer.writeIDref(ParameterParser.PARAMETER, operator.getName());
    }

    private void writeScaleOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(
                ScaleOperatorParser.SCALE_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>(ScaleOperatorParser.SCALE_FACTOR, operator.getTuning()),
                        getWeightAttribute(operator.getWeight())
                });
        writeParameter1Ref(writer, operator);
//        writeOperatorRef(writer, operator);
        writer.writeCloseTag(ScaleOperatorParser.SCALE_OPERATOR);
    }

    private void writeScaleOperator(Operator operator, XMLWriter writer, boolean indepedently) {
        writer.writeOpenTag(
                ScaleOperatorParser.SCALE_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>(ScaleOperatorParser.SCALE_FACTOR, operator.getTuning()),
                        getWeightAttribute(operator.getWeight()),
                        new Attribute.Default<String>(ScaleOperatorParser.SCALE_ALL_IND, indepedently ? "true" : "false")
                });
        writeParameter1Ref(writer, operator);
//        writeOperatorRef(writer, operator);
        writer.writeCloseTag(ScaleOperatorParser.SCALE_OPERATOR);
    }

    private void writeRandomWalkOperator(Operator operator, XMLWriter writer) {
        final String name = RandomWalkOperatorParser.RANDOM_WALK_OPERATOR;
        writer.writeOpenTag(
                name,
                new Attribute[]{
                        new Attribute.Default<Double>("windowSize", operator.getTuning()),
                        getWeightAttribute(operator.getWeight())
                });
        writeParameter1Ref(writer, operator);
//        writeOperatorRef(writer, operator);
        writer.writeCloseTag(name);
    }

    private void writeRandomWalkOperator(Operator operator, RandomWalkOperator.BoundaryCondition boundaryCondition, XMLWriter writer) {
        final String name = RandomWalkOperatorParser.RANDOM_WALK_OPERATOR;
        writer.writeOpenTag(
                name,
                new Attribute[]{
                        new Attribute.Default<Double>("windowSize", operator.getTuning()),
                        getWeightAttribute(operator.getWeight()),
                        new Attribute.Default<String>("boundaryCondition",
                                boundaryCondition.name()),
                        (operator.isAutoOptimize() == false ?
                                new Attribute.Default<Boolean>("autoOptimize", false) :
                                null)
                });
        writeParameter1Ref(writer, operator);
//        writeOperatorRef(writer, operator);
        writer.writeCloseTag(name);
    }

    private void writeIntegerRandomWalkOperator(Operator operator, XMLWriter writer) {

        int windowSize = (int) Math.round(operator.getTuning());
        if (windowSize < 1) windowSize = 1;
        final String name = RandomWalkIntegerOperatorParser.RANDOM_WALK_INTEGER_OPERATOR;
        writer.writeOpenTag(
                name,
                new Attribute[]{
                        new Attribute.Default<Integer>("windowSize", windowSize),
                        getWeightAttribute(operator.getWeight())
                });
        writeParameter1Ref(writer, operator);
//        writeOperatorRef(writer, operator);
        writer.writeCloseTag(name);
    }

    private void writeScaleAllOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(
                ScaleOperatorParser.SCALE_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>(ScaleOperatorParser.SCALE_FACTOR, operator.getTuning()),
                        new Attribute.Default<String>(ScaleOperatorParser.SCALE_ALL, "true"),
                        new Attribute.Default<String>(ScaleOperatorParser.IGNORE_BOUNDS, "true"),
                        getWeightAttribute(operator.getWeight())
                });

        if (operator.getParameter2() == null) {
            writeParameter1Ref(writer, operator);
        } else {
            writer.writeOpenTag(CompoundParameterParser.COMPOUND_PARAMETER);
            writeParameter1Ref(writer, operator);
//            writer.writeIDref(ParameterParser.PARAMETER, operator.getParameter2().getName());
            writeParameter2Ref(writer, operator);
            writer.writeCloseTag(CompoundParameterParser.COMPOUND_PARAMETER);
        }

        writer.writeCloseTag(ScaleOperatorParser.SCALE_OPERATOR);
    }

    private void writeDeltaOperator(Operator operator, boolean weighted, XMLWriter writer) {

        int[] parameterWeights = operator.getParameter1().getParameterDimensionWeights();
        Attribute[] attributes;

        if (weighted && parameterWeights != null && parameterWeights.length > 1) {
            String pw = "" + parameterWeights[0];
            for (int i = 1; i < parameterWeights.length; i++) {
                pw += " " + parameterWeights[i];
            }
            attributes = new Attribute[]{
                    new Attribute.Default<Double>(DeltaExchangeOperatorParser.DELTA, operator.getTuning()),
                    new Attribute.Default<String>(DeltaExchangeOperatorParser.PARAMETER_WEIGHTS, pw),
                    getWeightAttribute(operator.getWeight())
            };
        } else {
            attributes = new Attribute[]{
                    new Attribute.Default<Double>(DeltaExchangeOperatorParser.DELTA, operator.getTuning()),
                    getWeightAttribute(operator.getWeight())
            };
        }

        writer.writeOpenTag(DeltaExchangeOperatorParser.DELTA_EXCHANGE, attributes);
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(DeltaExchangeOperatorParser.DELTA_EXCHANGE);
    }

    private void writeIntegerDeltaOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(DeltaExchangeOperatorParser.DELTA_EXCHANGE,
                new Attribute[]{
                        new Attribute.Default<String>(DeltaExchangeOperatorParser.DELTA, Integer.toString((int) operator.getTuning())),
                        new Attribute.Default<String>("integer", "true"),
                        getWeightAttribute(operator.getWeight()),
                        new Attribute.Default<String>("autoOptimize", "false")
                }
        );
        writeParameter1Ref(writer, operator);
//        writeOperatorRef(writer, operator);
        writer.writeCloseTag(DeltaExchangeOperatorParser.DELTA_EXCHANGE);
    }

    private void writeSwapOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(SwapOperatorParser.SWAP_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<String>("size", Integer.toString((int) operator.getTuning())),
                        getWeightAttribute(operator.getWeight()),
                        new Attribute.Default<String>("autoOptimize", "false")
                }
        );
        writeParameter1Ref(writer, operator);
//        writeOperatorRef(writer, operator);
        writer.writeCloseTag(SwapOperatorParser.SWAP_OPERATOR);
    }

    private void writeBitFlipOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(BitFlipOperatorParser.BIT_FLIP_OPERATOR,
                getWeightAttribute(operator.getWeight()));
        writeParameter1Ref(writer, operator);
//        writeOperatorRef(writer, operator);
        writer.writeCloseTag(BitFlipOperatorParser.BIT_FLIP_OPERATOR);
    }

    private void writeBitFlipInSubstOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(BitFlipInSubstitutionModelOperator.BIT_FLIP_OPERATOR, new Attribute[]{
                new Attribute.Default<Double>(ScaleOperatorParser.SCALE_FACTOR, operator.getTuning()),
                getWeightAttribute(operator.getWeight())});
        writeParameter1Ref(writer, operator);
//        writeOperatorRef(writer, operator);
        PartitionSubstitutionModel model = (PartitionSubstitutionModel) operator.getOptions();
        writer.writeIDref(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL, model.getPrefix() + AbstractSubstitutionModel.MODEL);
        // <svsGeneralSubstitutionModel idref="originModel"/>
        writer.writeCloseTag(BitFlipInSubstitutionModelOperator.BIT_FLIP_OPERATOR);
    }

    private void writeRateBitExchangeOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(RateBitExchangeOperator.OPERATOR_NAME,
                getWeightAttribute(operator.getWeight()));

        writer.writeOpenTag(RateBitExchangeOperator.BITS);
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(RateBitExchangeOperator.BITS);

        writer.writeOpenTag(RateBitExchangeOperator.RATES);
        writeParameter2Ref(writer, operator);
        writer.writeCloseTag(RateBitExchangeOperator.RATES);

        writer.writeCloseTag(RateBitExchangeOperator.OPERATOR_NAME);
    }

    private void writeTreeBitMoveOperator(Operator operator, String treeModelPrefix, XMLWriter writer) {
        writer.writeOpenTag(TreeBitMoveOperatorParser.BIT_MOVE_OPERATOR,
                getWeightAttribute(operator.getWeight()));
        writer.writeIDref(DefaultTreeModel.TREE_MODEL, treeModelPrefix + DefaultTreeModel.TREE_MODEL);
        writer.writeCloseTag(TreeBitMoveOperatorParser.BIT_MOVE_OPERATOR);
    }

    private void writeUniformOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag("uniformOperator",
                getWeightAttribute(operator.getWeight()));
        writeParameter1Ref(writer, operator);
//        writeOperatorRef(writer, operator);
        writer.writeCloseTag("uniformOperator");
    }

    private void writeIntegerUniformOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(UniformIntegerOperatorParser.UNIFORM_INTEGER_OPERATOR,
                getWeightAttribute(operator.getWeight()));
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(UniformIntegerOperatorParser.UNIFORM_INTEGER_OPERATOR);
    }

    private void writeNarrowExchangeOperator(Operator operator, String treeModelPrefix, XMLWriter writer) {
        writer.writeOpenTag(ExchangeOperatorParser.NARROW_EXCHANGE,
                getWeightAttribute(operator.getWeight()));
        writer.writeIDref(DefaultTreeModel.TREE_MODEL, treeModelPrefix + DefaultTreeModel.TREE_MODEL);
        writer.writeCloseTag(ExchangeOperatorParser.NARROW_EXCHANGE);
    }

    private void writeWideExchangeOperator(Operator operator, String treeModelPrefix, XMLWriter writer) {
        writer.writeOpenTag(ExchangeOperatorParser.WIDE_EXCHANGE,
                getWeightAttribute(operator.getWeight()));
        writer.writeIDref(DefaultTreeModel.TREE_MODEL, treeModelPrefix + DefaultTreeModel.TREE_MODEL);
        writer.writeCloseTag(ExchangeOperatorParser.WIDE_EXCHANGE);
    }

    private void writeWilsonBaldingOperator(Operator operator, String treeModelPrefix, XMLWriter writer) {
        writer.writeOpenTag(WilsonBaldingParser.WILSON_BALDING,
                getWeightAttribute(operator.getWeight()));
        writer.writeIDref(DefaultTreeModel.TREE_MODEL, treeModelPrefix + DefaultTreeModel.TREE_MODEL);
        writer.writeCloseTag(WilsonBaldingParser.WILSON_BALDING);
    }

    private void writeSampleNonActiveOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(SampleNonActiveGibbsOperatorParser.SAMPLE_NONACTIVE_GIBBS_OPERATOR,
                getWeightAttribute(operator.getWeight()));

        writer.writeOpenTag(SampleNonActiveGibbsOperatorParser.DISTRIBUTION);
        writeOperatorRef(writer, operator);
        writer.writeCloseTag(SampleNonActiveGibbsOperatorParser.DISTRIBUTION);

        writer.writeOpenTag(SampleNonActiveGibbsOperatorParser.DATA_PARAMETER);
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(SampleNonActiveGibbsOperatorParser.DATA_PARAMETER);

        writer.writeOpenTag(SampleNonActiveGibbsOperatorParser.INDICATOR_PARAMETER);
        writeParameter2Ref(writer, operator);
        writer.writeCloseTag(SampleNonActiveGibbsOperatorParser.INDICATOR_PARAMETER);

        writer.writeCloseTag(SampleNonActiveGibbsOperatorParser.SAMPLE_NONACTIVE_GIBBS_OPERATOR);
    }

    private void writeSkyGridBlockUpdateOperator(Operator operator, String treePriorPrefix, XMLWriter writer) {
        writer.writeOpenTag(
                GMRFSkyrideBlockUpdateOperatorParser.GRID_BLOCK_UPDATE_OPERATOR,
                new Attribute[] {
                        new Attribute.Default<Double>(GMRFSkyrideBlockUpdateOperatorParser.SCALE_FACTOR, operator.getTuning()),
                        getWeightAttribute(operator.getWeight())
                }
        );
        writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD, treePriorPrefix + "skygrid");
        writer.writeCloseTag(GMRFSkyrideBlockUpdateOperatorParser.GRID_BLOCK_UPDATE_OPERATOR);
    }

    private void writeSkyGridHMCOperator(Operator operator, String treePriorPrefix, XMLWriter writer) {
        int nSteps = 50;
        double stepSize = 1E-2;
        String preconditioning = "none";
        int gradientCheckCount = 0;
        double gradientCheckTolerance = 1E-1;
        int preconditioningUpdateFrequency = 100;

        writer.writeOpenTag(
                HamiltonianMonteCarloOperatorParser.HMC_OPERATOR,
                new Attribute[]{
                        getWeightAttribute(operator.getWeight()),
                        new Attribute.Default<Integer>(HamiltonianMonteCarloOperatorParser.N_STEPS, nSteps),
                        new Attribute.Default<Double>(HamiltonianMonteCarloOperatorParser.STEP_SIZE, stepSize),
                        new Attribute.Default<String>(HamiltonianMonteCarloOperatorParser.MODE, "vanilla"),
                        new Attribute.Default<String>(AdaptableMCMCOperator.AUTO_OPTIMIZE, "true"),
                        new Attribute.Default<Integer>(HamiltonianMonteCarloOperatorParser.GRADIENT_CHECK_COUNT, gradientCheckCount),
                        new Attribute.Default<Double>(HamiltonianMonteCarloOperatorParser.GRADIENT_CHECK_TOLERANCE, gradientCheckTolerance),
                        new Attribute.Default<String>(HamiltonianMonteCarloOperatorParser.PRECONDITIONING, preconditioning),
                        new Attribute.Default<Integer>(HamiltonianMonteCarloOperatorParser.PRECONDITIONING_UPDATE_FREQUENCY, preconditioningUpdateFrequency)
                }
        );
        writer.writeIDref(CompoundGradientParser.COMPOUND_GRADIENT, treePriorPrefix + "full.skygrid.gradient");
        writer.writeIDref(CompoundParameterParser.COMPOUND_PARAMETER, treePriorPrefix + "skygrid.parameters");
        writer.writeOpenTag(
                SignTransformParser.NAME,
                new Attribute[]{
                        new Attribute.Default<Integer>(TransformParsers.START, 1),
                        new Attribute.Default<Integer>(TransformParsers.END, 1)
                }
        );
        writer.writeIDref(CompoundParameterParser.COMPOUND_PARAMETER, treePriorPrefix + "skygrid.parameters");
        writer.writeCloseTag(SignTransformParser.NAME);
        writer.writeCloseTag(HamiltonianMonteCarloOperatorParser.HMC_OPERATOR);
    }

    private void writeRelaxedClockHMCRateOperator(Operator operator, String prefix, XMLWriter writer) {
        int nSteps = 4;
        double stepSize = 1E-2;
        String preconditioning = "diagonal";
        int gradientCheckCount = 0;
        int preconditioningUpdateFrequency = 10;

        writer.writeOpenTag(
                HamiltonianMonteCarloOperatorParser.HMC_OPERATOR,
                new Attribute[]{
                        getWeightAttribute(operator.getWeight()),
                        new Attribute.Default<>(HamiltonianMonteCarloOperatorParser.N_STEPS, nSteps),
                        new Attribute.Default<>(HamiltonianMonteCarloOperatorParser.STEP_SIZE, stepSize),
                        new Attribute.Default<>("autoOptimize", true),
                        new Attribute.Default<>(HamiltonianMonteCarloOperatorParser.MODE, "vanilla"),
                        new Attribute.Default<>(HamiltonianMonteCarloOperatorParser.GRADIENT_CHECK_COUNT, gradientCheckCount),
                        new Attribute.Default<>(HamiltonianMonteCarloOperatorParser.PRECONDITIONING, preconditioning),
                        new Attribute.Default<>(HamiltonianMonteCarloOperatorParser.PRECONDITIONING_UPDATE_FREQUENCY, preconditioningUpdateFrequency)
                }
        );
        writer.writeOpenTag(JointGradientParser.JOINT_GRADIENT);
        writer.writeOpenTag(HessianWrapperParser.NAME);
        writer.writeIDref(DistributionLikelihood.DISTRIBUTION_LIKELIHOOD, prefix + BranchSpecificFixedEffects.RATES_PRIOR);
        writer.writeIDref(ParameterParser.PARAMETER, prefix + ClockType.HMC_CLOCK_BRANCH_RATES);
        writer.writeCloseTag(HessianWrapperParser.NAME);

        writer.writeOpenTag(BranchRateGradientParser.NAME, new Attribute.Default<>("traitName", "Sequence"));
        writer.writeIDref(TreeDataLikelihoodParser.TREE_DATA_LIKELIHOOD, prefix + "treeLikelihood");
        writer.writeCloseTag(BranchRateGradientParser.NAME);

        writer.writeCloseTag(JointGradientParser.JOINT_GRADIENT);

        writer.writeIDref(ParameterParser.PARAMETER, prefix + ClockType.HMC_CLOCK_BRANCH_RATES);

        writer.writeOpenTag(SignTransformParser.NAME);
        writer.writeIDref(ParameterParser.PARAMETER, prefix + ClockType.HMC_CLOCK_BRANCH_RATES);
        writer.writeCloseTag(SignTransformParser.NAME);
        writer.writeCloseTag(HamiltonianMonteCarloOperatorParser.HMC_OPERATOR);
    }

    private void writeRelaxedClockHMCScaleOperator(Operator operator, String prefix, XMLWriter writer) {
        int nSteps = 4;
        double stepSize = 1E-2;
        String preconditioning = "diagonal";
        int preconditioningUpdateFrequency = 10;
        int preconditioningDelay = 0;
        double drawVariance = 1.0;

        writer.writeOpenTag(
                HamiltonianMonteCarloOperatorParser.HMC_OPERATOR,
                new Attribute[]{
                        getWeightAttribute(operator.getWeight()),
                        new Attribute.Default<>(HamiltonianMonteCarloOperatorParser.N_STEPS, nSteps),
                        new Attribute.Default<>(HamiltonianMonteCarloOperatorParser.STEP_SIZE, stepSize),
                        new Attribute.Default<>("autoOptimize", true),
                        new Attribute.Default<>(HamiltonianMonteCarloOperatorParser.MODE, "vanilla"),
                        new Attribute.Default<>(HamiltonianMonteCarloOperatorParser.DRAW_VARIANCE, drawVariance),
                        new Attribute.Default<>(HamiltonianMonteCarloOperatorParser.PRECONDITIONING, preconditioning),
                        new Attribute.Default<>(HamiltonianMonteCarloOperatorParser.PRECONDITIONING_DELAY, preconditioningDelay),
                        new Attribute.Default<>(HamiltonianMonteCarloOperatorParser.PRECONDITIONING_UPDATE_FREQUENCY, preconditioningUpdateFrequency)
                }
        );
        writer.writeIDref(JointGradientParser.JOINT_GRADIENT, prefix + LocationScaleGradientParser.LOCATION_SCALE_JOINT_GRADIENT);
        writer.writeIDref(ParameterParser.PARAMETER, prefix + LocationScaleGradientParser.LOCATION_SCALE);
        writer.writeOpenTag(SignTransformParser.NAME);
        writer.writeIDref(ParameterParser.PARAMETER, prefix + LocationScaleGradientParser.LOCATION_SCALE);
        writer.writeCloseTag(SignTransformParser.NAME);
        writer.writeCloseTag(HamiltonianMonteCarloOperatorParser.HMC_OPERATOR);
    }

    private void writeShrinkageClockGibbsOperator(Operator operator, String prefix, XMLWriter writer) {
        writer.writeOpenTag(
                BayesianBridgeShrinkageOperatorParser.BAYESIAN_BRIDGE_PARSER,
                getWeightAttribute(operator.getWeight()));
        writer.writeIDref(AutoCorrelatedBranchRatesDistributionParser.AUTO_CORRELATED_RATES, prefix + "substBranchRatesPrior");
        writer.writeIDref(GAMMA_PRIOR, prefix + "globalScalePrior");
        writer.writeCloseTag(BayesianBridgeShrinkageOperatorParser.BAYESIAN_BRIDGE_PARSER);
    }

    private void writeShrinkageClockHMCOperator(Operator operator, String prefix, XMLWriter writer) {
        int nSteps = 5;
        double stepSize = 1E-2;
        int preconditioningUpdateFrequency = 1;
        double drawVariance = 1.0;
        int preconditioningDelay = 0;

        writer.writeOpenTag(
                HamiltonianMonteCarloOperatorParser.HMC_OPERATOR,
                new Attribute[]{
                        getWeightAttribute(operator.getWeight()),
                        new Attribute.Default<>(HamiltonianMonteCarloOperatorParser.N_STEPS, nSteps),
                        new Attribute.Default<>(HamiltonianMonteCarloOperatorParser.STEP_SIZE, stepSize),
                        new Attribute.Default<>(HamiltonianMonteCarloOperatorParser.MODE, "vanilla"),
                        new Attribute.Default<>(HamiltonianMonteCarloOperatorParser.DRAW_VARIANCE, drawVariance),
                        new Attribute.Default<>(HamiltonianMonteCarloOperatorParser.PRECONDITIONING_UPDATE_FREQUENCY, preconditioningUpdateFrequency),
                        new Attribute.Default<>(HamiltonianMonteCarloOperatorParser.PRECONDITIONING_DELAY, preconditioningDelay)
                }
        );
        writer.writeOpenTag(JointGradientParser.JOINT_GRADIENT);

        writer.writeComment("gradient of likelihood wrt increments");
        writer.writeOpenTag(BranchRateGradientWrtIncrementsParser.GRADIENT);

        writer.writeComment("gradient of likelihood wrt subst branch rates");
        writer.writeOpenTag(BranchRateGradientParser.NAME, new Attribute.Default<>("traitName", "Sequence"));
        writer.writeIDref(TreeDataLikelihoodParser.TREE_DATA_LIKELIHOOD, prefix + "treeLikelihood");
        writer.writeCloseTag(BranchRateGradientParser.NAME);

        writer.writeIDref(AutoCorrelatedGradientWrtIncrementsParser.GRADIENT, prefix + "incrementGradient");
        writer.writeCloseTag(BranchRateGradientWrtIncrementsParser.GRADIENT);

        writer.writeIDref(AutoCorrelatedGradientWrtIncrementsParser.GRADIENT, prefix + "incrementGradient");
        writer.writeCloseTag(JointGradientParser.JOINT_GRADIENT);

        writer.writeOpenTag("preconditioner");
        writer.writeIDref(BayesianBridgeDistributionModelParser.BAYESIAN_BRIDGE_DISTRIBUTION, prefix + "bbDistribution");
        writer.writeCloseTag("preconditioner");

        writer.writeCloseTag(HamiltonianMonteCarloOperatorParser.HMC_OPERATOR);
    }

    private void writeGMRFBlockUpdateOperator(Operator operator, String treePriorPrefix, XMLWriter writer) {
        writer.writeOpenTag(
                GMRFSkyrideBlockUpdateOperatorParser.BLOCK_UPDATE_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>(GMRFSkyrideBlockUpdateOperatorParser.SCALE_FACTOR, operator.getTuning()),
                        getWeightAttribute(operator.getWeight())
                }
        );
        writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD, treePriorPrefix + "skyride");
        writer.writeCloseTag(GMRFSkyrideBlockUpdateOperatorParser.BLOCK_UPDATE_OPERATOR);
    }

    private void writeScaleWithIndicatorsOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(
                ScaleOperatorParser.SCALE_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>(ScaleOperatorParser.SCALE_FACTOR, operator.getTuning()),
                        getWeightAttribute(operator.getWeight())
                });
        writeParameter1Ref(writer, operator);
        writer.writeOpenTag(ScaleOperatorParser.INDICATORS, new Attribute.Default<String>(ScaleOperatorParser.PICKONEPROB, "1.0"));
        writeParameter2Ref(writer, operator);
        writer.writeCloseTag(ScaleOperatorParser.INDICATORS);
        writer.writeCloseTag(ScaleOperatorParser.SCALE_OPERATOR);
    }

    // write multivariate operator

    private void writeSubtreeLeapOperator(Operator operator, String treeModelPrefix, XMLWriter writer) {
        writer.writeOpenTag(SubtreeLeapOperatorParser.SUBTREE_LEAP,
                new Attribute[]{
                        new Attribute.Default<Double>("size", operator.getTuning()),
                        getWeightAttribute(operator.getWeight())
                }
        );
        writer.writeIDref(DefaultTreeModel.TREE_MODEL, treeModelPrefix + DefaultTreeModel.TREE_MODEL);
        writer.writeCloseTag(SubtreeLeapOperatorParser.SUBTREE_LEAP);
    }

    private void writeFHSPROperator(Operator operator, String treeModelPrefix, XMLWriter writer) {
        writer.writeOpenTag(FixedHeightSubtreePruneRegraftOperatorParser.FIXED_HEIGHT_SUBTREE_PRUNE_REGRAFT,
                new Attribute[]{
                        getWeightAttribute(operator.getWeight())
                }
        );
        writer.writeIDref(DefaultTreeModel.TREE_MODEL, treeModelPrefix + DefaultTreeModel.TREE_MODEL);
        writer.writeCloseTag(FixedHeightSubtreePruneRegraftOperatorParser.FIXED_HEIGHT_SUBTREE_PRUNE_REGRAFT);
    }

    private void writeEmpiricalTreeSwapOperator(Operator operator, String treeModelPrefix, XMLWriter writer) {
        writer.writeOpenTag(EmpiricalTreeDistributionOperator.EMPIRICAL_TREE_DISTRIBUTION_OPERATOR,
                new Attribute[]{
                        getWeightAttribute(operator.getWeight())
                }
        );
        writer.writeIDref(EmpiricalTreeDistributionModel.EMPIRICAL_TREE_DISTRIBUTION_MODEL, treeModelPrefix + DefaultTreeModel.TREE_MODEL);
        writer.writeCloseTag(EmpiricalTreeDistributionOperator.EMPIRICAL_TREE_DISTRIBUTION_OPERATOR);
    }

    private void writeBFTUniformNodeHeightOperator(Operator operator, String treeModelPrefix, XMLWriter writer) {
        writer.writeOpenTag(NodeHeightOperatorParser.NODE_HEIGHT_OPERATOR, // should this be in the class and not the parser?
                new Attribute[]{
                    new Attribute.Default<String>("type",NodeHeightOperatorParser.OperatorType.UNIFORM.toString()),
                        getWeightAttribute(operator.getWeight())
                }
        );
        writer.writeIDref(ConstrainedTreeModel.CONSTRAINED_TREE_MODEL, treeModelPrefix + DefaultTreeModel.TREE_MODEL);
        writer.writeCloseTag(NodeHeightOperatorParser.NODE_HEIGHT_OPERATOR);
    }

    private void writeBFTRootScaleOperator(Operator operator, String treeModelPrefix, XMLWriter writer) {
        writer.writeOpenTag(NodeHeightOperatorParser.NODE_HEIGHT_OPERATOR, // should this be in the class and not the parser?
                new Attribute[]{
                    new Attribute.Default<String>("type",NodeHeightOperatorParser.OperatorType.SCALEROOT.toString()),
                    // todo scale factor here
                    new Attribute.Default<Double>(NodeHeightOperatorParser.SCALE_FACTOR, operator.getTuning()),
                        getWeightAttribute(operator.getWeight())
                }
        );
        writer.writeIDref(ConstrainedTreeModel.CONSTRAINED_TREE_MODEL, treeModelPrefix + DefaultTreeModel.TREE_MODEL);
        writer.writeCloseTag(NodeHeightOperatorParser.NODE_HEIGHT_OPERATOR);
    }

    private void writeUniformSPROperator(Operator operator, String treeModelPrefix, XMLWriter writer) {
        writer.writeOpenTag(UniformSubtreePruneRegraftParser.UNIFORM_SUBTREE_PRUNE_REGRAFT, // should this be in the class and not the parser?
            new Attribute[]{
                getWeightAttribute(operator.getWeight())
            });
        writer.writeIDref(ConstrainedTreeModel.CONSTRAINED_TREE_MODEL, treeModelPrefix + DefaultTreeModel.TREE_MODEL);
        writer.writeCloseTag(UniformSubtreePruneRegraftParser.UNIFORM_SUBTREE_PRUNE_REGRAFT);
    }


    // tuneable version of FHSPR but not currently being used
    private void writeSubtreeJumpOperator(Operator operator, String treeModelPrefix, XMLWriter writer) {
        writer.writeOpenTag(SubtreeJumpOperatorParser.SUBTREE_JUMP,
                new Attribute[]{
                        new Attribute.Default<Double>("size", operator.getTuning()),
                        getWeightAttribute(operator.getWeight())
                }
        );
        writer.writeIDref(DefaultTreeModel.TREE_MODEL, treeModelPrefix + DefaultTreeModel.TREE_MODEL);
        writer.writeCloseTag(SubtreeJumpOperatorParser.SUBTREE_JUMP);
    }

    private void writeSubtreeSlideOperator(Operator operator, String treeModelPrefix, XMLWriter writer) {
        writer.writeOpenTag(SubtreeSlideOperatorParser.SUBTREE_SLIDE,
                new Attribute[]{
                        new Attribute.Default<Double>("size", operator.getTuning()),
                        new Attribute.Default<String>("gaussian", "true"),
                        getWeightAttribute(operator.getWeight())
                }
        );
        writer.writeIDref(DefaultTreeModel.TREE_MODEL, treeModelPrefix + DefaultTreeModel.TREE_MODEL);
        writer.writeCloseTag(SubtreeSlideOperatorParser.SUBTREE_SLIDE);
    }

    private void writeUpDownOperator(String opTag, Operator operator, XMLWriter writer) {
        writer.writeOpenTag(opTag,
                new Attribute[]{
                        new Attribute.Default<Double>(ScaleOperatorParser.SCALE_FACTOR, operator.getTuning()),
                        getWeightAttribute(operator.getWeight())
                }
        );

        writer.writeOpenTag(UpDownOperatorParser.UP);
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(UpDownOperatorParser.UP);

        writer.writeOpenTag(UpDownOperatorParser.DOWN);
        writeParameter2Ref(writer, operator);
        writer.writeCloseTag(UpDownOperatorParser.DOWN);

        writer.writeCloseTag(opTag);
    }

    private void writeAdaptiveMultivariateOperator(Operator operator, XMLWriter writer) {

        if (NEW_AVMVN) {

            try {
                //determine how many parameters will be part of the AVMVN transition kernel
                int parameterCount = 0;
                for (Parameter parameter : options.selectParameters()) {
                    if (!parameter.isFixed() && parameter.isAdaptiveMultivariateCompatible) {
                        parameterCount++;
                    }
                }

                if (parameterCount > 0) {
                    //options set according to recommendations in AVMVN paper
                    int initial = 200 * parameterCount;
                    int burnin = initial / 2;

                    writer.writeOpenTag(AdaptableVarianceMultivariateNormalOperator.AVMVN_OPERATOR,
                            new Attribute[]{
                                    getWeightAttribute(operator.getWeight()),
                                    new Attribute.Default<Double>(AdaptableVarianceMultivariateNormalOperator.SCALE_FACTOR, operator.getTuning()),
                                    new Attribute.Default<Integer>(AdaptableVarianceMultivariateNormalOperator.INITIAL, initial),
                                    new Attribute.Default<Integer>(AdaptableVarianceMultivariateNormalOperator.BURNIN, burnin),
                                    new Attribute.Default<Double>(AdaptableVarianceMultivariateNormalOperator.BETA, 0.05),
                                    new Attribute.Default<Double>(AdaptableVarianceMultivariateNormalOperator.COEFFICIENT, 1.0),
                                    new Attribute.Default<Boolean>(AdaptableVarianceMultivariateNormalOperator.AUTO_OPTIMIZE, true),
                                    new Attribute.Default<Boolean>(AdaptableVarianceMultivariateNormalOperator.FORM_XTX, false)
                            });

                    ArrayList<Parameter> logList = new ArrayList<Parameter>();
                    ArrayList<Parameter> noList = new ArrayList<Parameter>();
                    ArrayList<Parameter> constrainedList = new ArrayList<Parameter>();

                    for (Parameter parameter : options.selectParameters()) {
                        if (parameter.isAdaptiveMultivariateCompatible && !parameter.isFixed()) {
                            //System.out.println(parameter.getName() + "   " + parameter.isMaintainedSum + " " + parameter.maintainedSum);
                            if (parameter.isNonNegative && !parameter.isMaintainedSum) {
                                logList.add(parameter);
                            } else if (parameter.isInRealSpace()) {
                                noList.add(parameter);
                            } else if (parameter.isMaintainedSum) {
                                constrainedList.add(parameter);
                            } else {
                                System.out.println("Parameter " + parameter + " should likely be equipped with a Dirichlet prior.");
                                System.out.println("Use of a Dirichlet prior for frequencies is set to: " + options.useNewFrequenciesPrior());
                                throw new UnsupportedOperationException("Parameter " + parameter.getName() + " with unidentified transformation.");
                            }
                        }
                    }

                    if (logList.size() > 0) {
                        writer.writeOpenTag(TransformParsers.TRANSFORM, new Attribute[]{new Attribute.Default<String>(TransformParsers.TYPE, new Transform.LogTransform().getTransformName())});
                        for (Parameter parameter : logList) {
                            writer.writeIDref(ParameterParser.PARAMETER, parameter.getName());
                        }
                        writer.writeCloseTag(TransformParsers.TRANSFORM);
                    }

                    if (noList.size() > 0) {
                        writer.writeOpenTag(TransformParsers.TRANSFORM, new Attribute[]{new Attribute.Default<String>(TransformParsers.TYPE, new Transform.NoTransform().getTransformName())});
                        for (Parameter parameter : noList) {
                            writer.writeIDref(ParameterParser.PARAMETER, parameter.getName());
                        }
                        writer.writeCloseTag(TransformParsers.TRANSFORM);
                    }

                    for (Parameter parameter : constrainedList) {
                        writer.writeOpenTag(TransformParsers.TRANSFORM, new Attribute[]{new Attribute.Default<String>(TransformParsers.TYPE, new Transform.LogConstrainedSumTransform().getTransformName()),
                                new Attribute.Default<Double>(TransformParsers.SUM, parameter.maintainedSum)});
                        writer.writeIDref(ParameterParser.PARAMETER, parameter.getName());
                        writer.writeCloseTag(TransformParsers.TRANSFORM);
                    }

                    writer.writeCloseTag(AdaptableVarianceMultivariateNormalOperator.AVMVN_OPERATOR);

                }

            } catch (UnsupportedOperationException unSup) {
                System.out.println(unSup);
            }

        } else {

            //determine how many parameters will be part of the AVMVN transition kernel
            int parameterCount = 0;
            for (Parameter parameter : options.selectParameters()) {
                if (!parameter.isFixed() && parameter.isAdaptiveMultivariateCompatible) {
                    parameterCount++;
                }
            }

            if (parameterCount > 0) {

                //options set according to recommendations in AVMVN paper
                int initial = 200 * parameterCount;
                int burnin = initial / 2;

                writer.writeOpenTag(AdaptableVarianceMultivariateNormalOperator.AVMVN_OPERATOR,
                        new Attribute[]{
                                getWeightAttribute(operator.getWeight()),
                                new Attribute.Default<Double>(AdaptableVarianceMultivariateNormalOperator.SCALE_FACTOR, operator.getTuning()),
                                new Attribute.Default<Integer>(AdaptableVarianceMultivariateNormalOperator.INITIAL, initial),
                                new Attribute.Default<Integer>(AdaptableVarianceMultivariateNormalOperator.BURNIN, burnin),
                                new Attribute.Default<Double>(AdaptableVarianceMultivariateNormalOperator.BETA, 0.05),
                                new Attribute.Default<Double>(AdaptableVarianceMultivariateNormalOperator.COEFFICIENT, 1.0),
                                new Attribute.Default<Boolean>(AdaptableVarianceMultivariateNormalOperator.AUTO_OPTIMIZE, true),
                                new Attribute.Default<Boolean>(AdaptableVarianceMultivariateNormalOperator.FORM_XTX, false)
                        });

                // @todo Need to collate only the parameters being controlled by this here.

                for (Parameter parameter : options.selectParameters()) {
                    if (parameter.isAdaptiveMultivariateCompatible) {
                        writer.writeIDref(ParameterParser.PARAMETER, parameter.getName());
                    }
                }

                //set appropriate transformations for all parameters
                //TODO: we should aggregate as best as possible the different transformation so as to have fewer attributes
                int startTransform = 0;
                for (Parameter parameter : options.selectParameters()) {
                    if (parameter.isAdaptiveMultivariateCompatible) {
                        if (parameter.isNonNegative) {
                            writer.writeTag(TransformParsers.TRANSFORM, new Attribute[]{new Attribute.Default<String>(TransformParsers.TYPE, new Transform.LogTransform().getTransformName()),
                                    new Attribute.Default<Integer>(TransformParsers.START, startTransform),
                                    new Attribute.Default<Integer>(TransformParsers.END, startTransform + parameter.getDimensionWeight()),
                            }, true);
                            startTransform += parameter.getDimensionWeight();
                            System.out.println(parameter + ": " + parameter.getDimensionWeight());

                        } else { // -Inf to Inf
                            writer.writeTag(TransformParsers.TRANSFORM, new Attribute[]{new Attribute.Default<String>(TransformParsers.TYPE, new Transform.NoTransform().getTransformName()),
                                    new Attribute.Default<Integer>(TransformParsers.START, startTransform),
                                    new Attribute.Default<Integer>(TransformParsers.END, startTransform + parameter.getDimensionWeight()),
                            }, true);
                            startTransform += parameter.getDimensionWeight();
                            System.out.println(parameter + ": " + parameter.getDimensionWeight());
                        }
                    }
                }

                writer.writeCloseTag(AdaptableVarianceMultivariateNormalOperator.AVMVN_OPERATOR);

            }

        }
    }

    private Attribute getWeightAttribute(double weight) {
        if (weight == (int) weight) {
            return new Attribute.Default<Integer>("weight", (int) weight);
        } else {
            return new Attribute.Default<Double>("weight", weight);
        }
    }
}
