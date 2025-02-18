/*
 * TreePriorGenerator.java
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
import dr.app.beauti.options.Parameter;
import dr.app.beauti.options.PartitionTreeModel;
import dr.app.beauti.options.PartitionTreePrior;
import dr.app.beauti.types.StartingTreeType;
import dr.app.beauti.types.TreePriorParameterizationType;
import dr.app.beauti.types.TreePriorType;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.util.Taxa;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.GMRFSkyrideGradient;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodelxml.coalescent.CoalescentLikelihoodParser;
import dr.evomodelxml.coalescent.GMRFSkyrideGradientParser;
import dr.evomodelxml.coalescent.GMRFSkyrideLikelihoodParser;
import dr.evomodelxml.coalescent.demographicmodel.ConstantPopulationModelParser;
import dr.evomodelxml.coalescent.demographicmodel.ExpansionModelParser;
import dr.evomodelxml.coalescent.demographicmodel.ExponentialGrowthModelParser;
import dr.evomodelxml.coalescent.demographicmodel.LogisticGrowthModelParser;
import dr.evomodelxml.speciation.BirthDeathModelParser;
import dr.evomodelxml.speciation.BirthDeathSerialSamplingModelParser;
import dr.evomodelxml.speciation.SpeciationLikelihoodParser;
import dr.evomodelxml.speciation.YuleModelParser;
import dr.evoxml.TaxaParser;
import dr.inference.model.ParameterParser;
import dr.inferencexml.distribution.GammaDistributionModelParser;
import dr.inferencexml.distribution.PriorParsers;
import dr.inferencexml.hmc.CompoundGradientParser;
import dr.inferencexml.hmc.GradientWrapperParser;
import dr.inferencexml.hmc.JointGradientParser;
import dr.inferencexml.model.CompoundParameterParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class TreePriorGenerator extends Generator {

    public TreePriorGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

//    void writeTreePrior(PartitionTreePrior prior, PartitionTreeModel model, XMLWriter writer) {    // for species, partitionName.treeModel
//        setModelPrefix(prior.getPrefix()); // only has prefix, if (options.getPartitionTreePriors().size() > 1)
//
//        writePriorLikelihood(prior, model, writer);
//    }

    /**
     * Write a tree prior (coalescent or speciational) model
     *
     * @param prior  the partition tree prior
     * @param writer the writer
     */
    void writeTreePriorModel(PartitionTreePrior prior, XMLWriter writer) {

        String prefix = prior.getPrefix();

        String initialPopSize = null;

        TreePriorType nodeHeightPrior = prior.getNodeHeightPrior();
        Units.Type units = options.units;
        TreePriorParameterizationType parameterization = prior.getParameterization();

        switch (nodeHeightPrior) {
            case CONSTANT:
                writer.writeComment("A prior assumption that the population size has remained constant",
                        "throughout the time spanned by the genealogy.");
                writer.writeOpenTag(
                        ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + "constant"),
                                new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(options.units))
                        }
                );

                writer.writeOpenTag(ConstantPopulationModelParser.POPULATION_SIZE);
                writeParameter("constant.popSize", prior, writer);
                writer.writeCloseTag(ConstantPopulationModelParser.POPULATION_SIZE);
                writer.writeCloseTag(ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL);

                break;

            case EXPONENTIAL:
                // generate an exponential prior tree

                writer.writeComment("A prior assumption that the population size has grown exponentially",
                        "throughout the time spanned by the genealogy.");
                writer.writeOpenTag(
                        ExponentialGrowthModelParser.EXPONENTIAL_GROWTH_MODEL,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + "exponential"),
                                new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(options.units))
                        }
                );

                // write pop size socket
                writer.writeOpenTag(ExponentialGrowthModelParser.POPULATION_SIZE);
                writeParameter("exponential.popSize", prior, writer);
                writer.writeCloseTag(ExponentialGrowthModelParser.POPULATION_SIZE);

                if (parameterization == TreePriorParameterizationType.GROWTH_RATE) {
                    // write growth rate socket
                    writer.writeOpenTag(ExponentialGrowthModelParser.GROWTH_RATE);
                    writeParameter("exponential.growthRate", prior, writer);
                    writer.writeCloseTag(ExponentialGrowthModelParser.GROWTH_RATE);
                } else {
                    // write doubling time socket
                    writer.writeOpenTag(ExponentialGrowthModelParser.DOUBLING_TIME);
                    writeParameter("exponential.doublingTime", prior, writer);
                    writer.writeCloseTag(ExponentialGrowthModelParser.DOUBLING_TIME);
                }

                writer.writeCloseTag(ExponentialGrowthModelParser.EXPONENTIAL_GROWTH_MODEL);

                break;

            case LOGISTIC:
                // generate an exponential prior tree
                writer.writeComment("A prior assumption that the population size has grown logistically",
                        "throughout the time spanned by the genealogy.");
                writer.writeOpenTag(
                        LogisticGrowthModelParser.LOGISTIC_GROWTH_MODEL,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + "logistic"),
                                new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(options.units))
                        }
                );

                // write pop size socket
                writer.writeOpenTag(LogisticGrowthModelParser.POPULATION_SIZE);
                writeParameter("logistic.popSize", prior, writer);
                writer.writeCloseTag(LogisticGrowthModelParser.POPULATION_SIZE);

                if (parameterization == TreePriorParameterizationType.GROWTH_RATE) {
                    // write growth rate socket
                    writer.writeOpenTag(LogisticGrowthModelParser.GROWTH_RATE);
                    writeParameter("logistic.growthRate", prior, writer);
                    writer.writeCloseTag(LogisticGrowthModelParser.GROWTH_RATE);
                } else {
                    // write doubling time socket
                    writer.writeOpenTag(LogisticGrowthModelParser.DOUBLING_TIME);
                    writeParameter("logistic.doublingTime", prior, writer);
                    writer.writeCloseTag(LogisticGrowthModelParser.DOUBLING_TIME);
                }

                // write logistic t50 socket
                writer.writeOpenTag(LogisticGrowthModelParser.TIME_50);

//	            if (options.clockModelOptions.getRateOptionClockModel() == FixRateType.FIX_MEAN
//	        			|| options.clockModelOptions.getRateOptionClockModel() == FixRateType.RELATIVE_TO) {
//            		writer.writeComment("No calibration");
//                writer.writeComment("logistic.t50 initial always has to < treeRootHeight initial");
//                dr.app.beauti.options.Parameter priorPara = prior.getParameter("logistic.t50");
//
//                double initRootHeight;
//                if (options.isShareSameTreePrior()) {
//                    initRootHeight = priorPara.initial;
//                    for (PartitionTreeModel tree : options.getPartitionTreeModels()) {
//                        double tmpRootHeight = tree.getParameter("treeModel.rootHeight").initial;
//                        if (initRootHeight > tmpRootHeight) { // take min
//                            initRootHeight = tmpRootHeight;
//                        }
//                    }
//                } else {
//                    initRootHeight = prior.getTreeModel().getParameter("treeModel.rootHeight").initial;
//                }
//                // logistic.t50 initial always has to < treeRootHeight initial
//                if (priorPara.initial >= initRootHeight) {
//                    priorPara.initial = initRootHeight / 2; // tree prior.initial has to < treeRootHeight.initial
//                }
//	            } else {
//	            	writer.writeComment("Has calibration");
//
//	            	throw new IllegalArgumentException("This function is not available in this release !");
//	            }

                writeParameter("logistic.t50", prior, writer);
                writer.writeCloseTag(LogisticGrowthModelParser.TIME_50);

                writer.writeCloseTag(LogisticGrowthModelParser.LOGISTIC_GROWTH_MODEL);

                initialPopSize = "logistic.popSize";

                break;

            case EXPANSION:
                // generate an exponential prior tree
                writer.writeComment("A prior assumption that the population size has grown exponentially",
                        "from some ancestral population size in the past.");
                writer.writeOpenTag(
                        ExpansionModelParser.EXPANSION_MODEL,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + "expansion"),
                                new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(options.units))
                        }
                );

                // write pop size socket
                writeParameter(ExpansionModelParser.POPULATION_SIZE, "expansion.popSize", prior, writer);

                if (parameterization == TreePriorParameterizationType.GROWTH_RATE) {
                    // write growth rate socket
                    writeParameter(ExpansionModelParser.GROWTH_RATE, "expansion.growthRate", prior, writer);
                } else {
                    // write doubling time socket
                    writeParameter(ExpansionModelParser.DOUBLING_TIME, "expansion.doublingTime", prior, writer);
                }

                // write ancestral proportion socket
                writeParameter(ExpansionModelParser.ANCESTRAL_POPULATION_PROPORTION, "expansion.ancestralProportion", prior, writer);

                writer.writeCloseTag(ExpansionModelParser.EXPANSION_MODEL);

                initialPopSize = "expansion.popSize";

                break;

            case YULE:
            case YULE_CALIBRATION:
                if (nodeHeightPrior == TreePriorType.YULE_CALIBRATION) {
                    writer.writeComment("Calibrated Yule: Heled J, Drummond AJ (2011), Syst Biol, doi: " +
                            "10.1093/sysbio/syr087");
                } else {
                    writer.writeComment("A prior on the distribution node heights defined given",
                            "a Yule speciation process (a pure birth process).");
                }
                writer.writeOpenTag(
                        YuleModelParser.YULE_MODEL,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + YuleModelParser.YULE),
                                new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
                        }
                );

                writeParameter(YuleModelParser.BIRTH_RATE, "yule.birthRate", prior, writer);
                writer.writeCloseTag(YuleModelParser.YULE_MODEL);

                break;

            case BIRTH_DEATH:
            case BIRTH_DEATH_INCOMPLETE_SAMPLING:
                writer.writeComment("A prior on the distribution node heights defined given");
                writer.writeComment(nodeHeightPrior == TreePriorType.BIRTH_DEATH_INCOMPLETE_SAMPLING ?
                        BirthDeathModelParser.getCitationRHO() : BirthDeathModelParser.getCitation());
                writer.writeOpenTag(
                        BirthDeathModelParser.BIRTH_DEATH_MODEL,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + BirthDeathModelParser.BIRTH_DEATH),
                                new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
                        }
                );

                writeParameter(BirthDeathModelParser.BIRTHDIFF_RATE, BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME, prior, writer);
                writeParameter(BirthDeathModelParser.RELATIVE_DEATH_RATE, BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME, prior, writer);

                if (nodeHeightPrior == TreePriorType.BIRTH_DEATH_INCOMPLETE_SAMPLING) {
                    writeParameter(BirthDeathModelParser.SAMPLE_PROB,
                            BirthDeathModelParser.BIRTH_DEATH + "." + BirthDeathModelParser.SAMPLE_PROB, prior, writer);
                }

                writer.writeCloseTag(BirthDeathModelParser.BIRTH_DEATH_MODEL);

                break;

            case BIRTH_DEATH_SERIAL_SAMPLING:
                writer.writeComment(BirthDeathSerialSamplingModelParser.getCitationPsiOrg());

                writer.writeOpenTag(
                        BirthDeathSerialSamplingModelParser.BIRTH_DEATH_SERIAL_MODEL,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + BirthDeathSerialSamplingModelParser.BDSS),
                                new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units)),
                                new Attribute.Default<Boolean>(BirthDeathSerialSamplingModelParser.HAS_FINAL_SAMPLE, false)
                        }
                );

                writeParameter(BirthDeathSerialSamplingModelParser.LAMBDA,
                        BirthDeathSerialSamplingModelParser.BDSS + "." + BirthDeathSerialSamplingModelParser.LAMBDA, prior, writer);
                writeParameter(BirthDeathSerialSamplingModelParser.RELATIVE_MU,
                        BirthDeathSerialSamplingModelParser.BDSS + "." + BirthDeathSerialSamplingModelParser.RELATIVE_MU, prior, writer);
//                writeParameter(BirthDeathSerialSamplingModelParser.SAMPLE_PROBABILITY,
//                        BirthDeathSerialSamplingModelParser.BDSS + "." + BirthDeathSerialSamplingModelParser.SAMPLE_PROBABILITY, prior, writer);
                writeParameter(BirthDeathSerialSamplingModelParser.PSI,
                        BirthDeathSerialSamplingModelParser.BDSS + "." + BirthDeathSerialSamplingModelParser.PSI, prior, writer);
                writeParameter(BirthDeathSerialSamplingModelParser.ORIGIN,
                        BirthDeathSerialSamplingModelParser.BDSS + "." + BirthDeathSerialSamplingModelParser.ORIGIN, prior, writer);

                writer.writeCloseTag(BirthDeathSerialSamplingModelParser.BIRTH_DEATH_SERIAL_MODEL);

                break;

//            case BIRTH_DEATH_BASIC_REPRODUCTIVE_NUMBER:
//                writer.writeComment(BirthDeathSerialSamplingModelParser.getCitationRT());
//
//                writer.writeOpenTag(
//                        BirthDeathEpidemiologyModelParser.BIRTH_DEATH_EPIDEMIOLOGY,
//                        new Attribute[]{
//                                new Attribute.Default<String>(XMLParser.ID, prefix + BirthDeathEpidemiologyModelParser.BIRTH_DEATH_EPIDEMIOLOGY),
//                                new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
//                        }
//                );
//
//                writeParameter(BirthDeathEpidemiologyModelParser.R0,
//                        BirthDeathEpidemiologyModelParser.R0, prior, writer);
//                writeParameter(BirthDeathEpidemiologyModelParser.RECOVERY_RATE,
//                        BirthDeathEpidemiologyModelParser.RECOVERY_RATE, prior, writer);
//                writeParameter(BirthDeathEpidemiologyModelParser.SAMPLING_PROBABILITY,
//                        BirthDeathEpidemiologyModelParser.SAMPLING_PROBABILITY, prior, writer);
//                writeParameter(BirthDeathEpidemiologyModelParser.ORIGIN,
//                        BirthDeathEpidemiologyModelParser.ORIGIN, prior, writer);
//
//                writer.writeCloseTag(BirthDeathEpidemiologyModelParser.BIRTH_DEATH_EPIDEMIOLOGY);
//
//                break;
            case SKYGRID_HMC:
            case SKYGRID:
            case GMRF_SKYRIDE:
                // do nothing here...
                break;
            default:
                throw new UnsupportedOperationException("Unknown Tree Prior type");
        }

        if (nodeHeightPrior != TreePriorType.CONSTANT && nodeHeightPrior != TreePriorType.EXPONENTIAL) {
            // If the node height prior is not one of these two then we need to simulate a
            // random starting tree under a constant size coalescent.
            writer.writeComment("This is a simple constant population size coalescent model",
                    "that is used to generate an initial tree for the chain.");
            writer.writeOpenTag(
                    ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, prefix + "initialDemo"),
                            new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
                    }
            );

            writer.writeOpenTag(ConstantPopulationModelParser.POPULATION_SIZE);
            if (initialPopSize != null) {
                writer.writeIDref(ParameterParser.PARAMETER, prefix + initialPopSize);
            } else {
                writeParameter(prefix + "initialDemo.popSize", 1, 100.0, Double.NaN, Double.NaN, writer);
            }
            writer.writeCloseTag(ConstantPopulationModelParser.POPULATION_SIZE);
            writer.writeCloseTag(ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL);
        }

//        if (nodeHeightPrior == TreePriorType.BIRTH_DEATH_BASIC_REPRODUCTIVE_NUMBER) {
//            writer.writeComment("R0 = b/(b*d+s*r)");
//            writer.writeOpenTag(RPNcalculatorStatisticParser.RPN_STATISTIC,
//                    new Attribute[]{
//                            new Attribute.Default<String>(XMLParser.ID, modelPrefix + "R0")
//                    });
//
//            writer.writeOpenTag(RPNcalculatorStatisticParser.VARIABLE,
//                    new Attribute[]{
//                            new Attribute.Default<String>(Statistic.NAME, modelPrefix + "b")
//                    });
//            writeParameterRef(modelPrefix + BirthDeathSerialSamplingModelParser.BDSS + "." + BirthDeathSerialSamplingModelParser.LAMBDA, writer);
//            writer.writeCloseTag(RPNcalculatorStatisticParser.VARIABLE);
//
//            writer.writeOpenTag(RPNcalculatorStatisticParser.VARIABLE,
//                    new Attribute[]{
//                            new Attribute.Default<String>(Statistic.NAME, modelPrefix + "d")
//                    });
//            writeParameterRef(modelPrefix + BirthDeathSerialSamplingModelParser.BDSS + "." + BirthDeathSerialSamplingModelParser.RELATIVE_MU, writer);
//            writer.writeCloseTag(RPNcalculatorStatisticParser.VARIABLE);
//
//            writer.writeOpenTag(RPNcalculatorStatisticParser.VARIABLE,
//                    new Attribute[]{
//                            new Attribute.Default<String>(Statistic.NAME, modelPrefix + "s")
//                    });
//            writeParameterRef(modelPrefix + BirthDeathSerialSamplingModelParser.BDSS + "." + BirthDeathSerialSamplingModelParser.PSI, writer);
//            writer.writeCloseTag(RPNcalculatorStatisticParser.VARIABLE);
//
//            writer.writeOpenTag(RPNcalculatorStatisticParser.VARIABLE,
//                    new Attribute[]{
//                            new Attribute.Default<String>(Statistic.NAME, modelPrefix + "r")
//                    });
//            writeParameterRef(modelPrefix + BirthDeathSerialSamplingModelParser.BDSS + "." + BirthDeathSerialSamplingModelParser.R, writer);
//            writer.writeCloseTag(RPNcalculatorStatisticParser.VARIABLE);
//
//            writer.writeOpenTag(RPNcalculatorStatisticParser.EXPRESSION,
//                    new Attribute[]{
//                            new Attribute.Default<String>(Statistic.NAME, modelPrefix + "R0")
//                    });
//            writer.writeText(modelPrefix + "b " + modelPrefix + "b " + modelPrefix + "d " + "* " + modelPrefix + "s " + modelPrefix + "r " + "* + /");
//            writer.writeCloseTag(RPNcalculatorStatisticParser.EXPRESSION);
//
//            writer.writeCloseTag(RPNcalculatorStatisticParser.RPN_STATISTIC);
//        }
    }

    /**
     * Write the prior on node heights (coalescent or speciational models)
     *
     * @param model  PartitionTreeModel
     * @param writer the writer
     */
    void writePriorLikelihood(PartitionTreeModel model, XMLWriter writer) {

        //tree model prefix
        String prefix = model.getPrefix();

        PartitionTreePrior prior = model.getPartitionTreePrior();
        TreePriorType treePrior = prior.getNodeHeightPrior();

//        String priorPrefix = prior.getPrefix();

        switch (treePrior) {
            case YULE:
            case BIRTH_DEATH:
            case BIRTH_DEATH_INCOMPLETE_SAMPLING:
            case BIRTH_DEATH_SERIAL_SAMPLING:
//            case BIRTH_DEATH_BASIC_REPRODUCTIVE_NUMBER:
            case YULE_CALIBRATION:
                // generate a speciational process
                writer.writeComment("Generate a speciation likelihood for Yule or Birth Death");
                writer.writeOpenTag(
                        SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + "speciation")
                        }
                );

                // write pop size socket
                writer.writeOpenTag(SpeciationLikelihoodParser.MODEL);
                writeNodeHeightPriorModelRef(prior, writer);
                writer.writeCloseTag(SpeciationLikelihoodParser.MODEL);
                writer.writeOpenTag(SpeciationLikelihoodParser.TREE);
                writer.writeIDref(DefaultTreeModel.TREE_MODEL, prefix + DefaultTreeModel.TREE_MODEL);
                writer.writeCloseTag(SpeciationLikelihoodParser.TREE);

                if (treePrior == TreePriorType.YULE_CALIBRATION) {

                    if (options.treeModelOptions.isNodeCalibrated(model) == 0) {
                        writer.writeOpenTag(SpeciationLikelihoodParser.CALIBRATION,
                                new Attribute[]{
                                        new Attribute.Default<String>(SpeciationLikelihoodParser.CORRECTION, prior.getCalibCorrectionType().toString())
                                });
                        writer.writeOpenTag(SpeciationLikelihoodParser.POINT);

                        String taxaId;
                        if (options.hasIdenticalTaxa()) {
                            taxaId = TaxaParser.TAXA;
                        } else {
                            taxaId = options.getDataPartitions(model).get(0).getPrefix() + TaxaParser.TAXA;
                        }
                        writer.writeIDref(TaxaParser.TAXA, taxaId);

                        writeDistribution(model.getParameter("treeModel.rootHeight"), true, writer);

                        writer.writeCloseTag(SpeciationLikelihoodParser.POINT);
                        writer.writeCloseTag(SpeciationLikelihoodParser.CALIBRATION);

                    } else if (options.treeModelOptions.isNodeCalibrated(model) == 1) {
                        // should be only 1 calibrated internal node with monophyletic for each tree at moment
                        Taxa t = (Taxa) options.getKeysFromValue(options.taxonSetsTreeModel, model).get(0);
                        Parameter nodeCalib = options.getStatistic(t);

                        writer.writeOpenTag(SpeciationLikelihoodParser.CALIBRATION,
                                new Attribute[]{
                                        new Attribute.Default<String>(SpeciationLikelihoodParser.CORRECTION, prior.getCalibCorrectionType().toString())
                                });
                        writer.writeOpenTag(SpeciationLikelihoodParser.POINT);

                        writer.writeIDref(TaxaParser.TAXA, t.getId());
                        writeDistribution(nodeCalib, true, writer);

                        writer.writeCloseTag(SpeciationLikelihoodParser.POINT);
                        writer.writeCloseTag(SpeciationLikelihoodParser.CALIBRATION);

                        if (!options.treeModelOptions.isNodeCalibrated(nodeCalib)) {
                            throw new IllegalArgumentException("Calibrated Yule model requires a calibration to be specified for node, " +
                                    nodeCalib.getName() + ".");
                        }
                    }
                }

                writer.writeCloseTag(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD);

                break;

            // AR - removed this special case for Logistic - it causes terrible problems.
            // Need to put informative priors on Logistic parameters.
//    		case LOGISTIC:
//    			writer.writeComment("Generate a boolean likelihood for Coalescent: Logistic Growth");
//    			writer.writeOpenTag(
//    	                BooleanLikelihoodParser.BOOLEAN_LIKELIHOOD,
//    	                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, modelPrefix + "booleanLikelihood1")}
//    	        );
//    	        writer.writeOpenTag(
//    	                TestStatisticParser.TEST_STATISTIC,
//    	                new Attribute[]{
//    	                        new Attribute.Default<String>(XMLParser.ID, "test1"),
//    	                        new Attribute.Default<String>(Statistic.NAME, "test1")
//    	                }
//    	        );
//    	        writer.writeIDref(ParameterParser.PARAMETER, priorPrefix + "logistic.t50");
//    	        writer.writeOpenTag("lessThan");
//    	        writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "treeModel.rootHeight");
//    	        writer.writeCloseTag("lessThan");
//    	        writer.writeCloseTag(TestStatisticParser.TEST_STATISTIC);
//    	        writer.writeCloseTag(BooleanLikelihoodParser.BOOLEAN_LIKELIHOOD);
//
//    	        writer.writeOpenTag(
//	                    CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD,
//	                    new Attribute[]{new Attribute.Default<String>(XMLParser.ID, modelPrefix + COALESCENT)}
//	            );
//	            writer.writeOpenTag(CoalescentLikelihoodParser.MODEL);
//	            writeNodeHeightPriorModelRef(prior, writer);
//	            writer.writeCloseTag(CoalescentLikelihoodParser.MODEL);
//	            writer.writeOpenTag(CoalescentLikelihoodParser.POPULATION_TREE);
//	            writer.writeIDref(DefaultTreeModel.TREE_MODEL, modelPrefix + DefaultTreeModel.TREE_MODEL);
//	            writer.writeCloseTag(CoalescentLikelihoodParser.POPULATION_TREE);
//	            writer.writeCloseTag(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD);
//
//    	        break;

            case GMRF_SKYRIDE:
                writer.writeComment("Generate a gmrfSkyrideLikelihood for GMRF Bayesian Skyride process");
                writer.writeOpenTag(
                        GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + "skyride"),
                                new Attribute.Default<String>(GMRFSkyrideLikelihoodParser.TIME_AWARE_SMOOTHING,
                                        prior.getSkyrideSmoothing() == TreePriorParameterizationType.TIME_AWARE_SKYRIDE ? "true" : "false"),
                                new Attribute.Default<String>(GMRFSkyrideLikelihoodParser.RANDOMIZE_TREE,
                                        //TODO For GMRF, tree model/tree prior combination not implemented by BEAST yet. The validation is in BeastGenerator.checkOptions()
                                        options.getPartitionTreeModels(prior).get(0).getStartingTreeType() == StartingTreeType.UPGMA ? "true" : "false"),
                        }
                );

                int skyrideIntervalCount = options.taxonList.getTaxonCount() - 1;
                writer.writeOpenTag(GMRFSkyrideLikelihoodParser.POPULATION_PARAMETER);
                writer.writeComment("skyride.logPopSize is in log units unlike other popSize");
                writeParameter(prior.getParameter("skyride.logPopSize"), skyrideIntervalCount, writer);
                writer.writeCloseTag(GMRFSkyrideLikelihoodParser.POPULATION_PARAMETER);

                writer.writeOpenTag(GMRFSkyrideLikelihoodParser.GROUP_SIZES);
                writeParameter(prior.getParameter("skyride.groupSize"), skyrideIntervalCount, writer);
                writer.writeCloseTag(GMRFSkyrideLikelihoodParser.GROUP_SIZES);

                writer.writeOpenTag(GMRFSkyrideLikelihoodParser.PRECISION_PARAMETER);
                writeParameter(prior.getParameter("skyride.precision"), 1, writer);
                writer.writeCloseTag(GMRFSkyrideLikelihoodParser.PRECISION_PARAMETER);

                writer.writeOpenTag(GMRFSkyrideLikelihoodParser.POPULATION_TREE);
                writer.writeIDref(DefaultTreeModel.TREE_MODEL, prefix + DefaultTreeModel.TREE_MODEL);
                writer.writeCloseTag(GMRFSkyrideLikelihoodParser.POPULATION_TREE);

                writer.writeCloseTag(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD);

                break;

            case SKYGRID:
            case SKYGRID_HMC:
                break;

            default:
                // generate a coalescent process
                writer.writeComment("Generate a coalescent likelihood");
                writer.writeOpenTag(
                        CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD,
                        new Attribute[]{new Attribute.Default<String>(XMLParser.ID, prefix + COALESCENT)}
                );
                writer.writeOpenTag(CoalescentLikelihoodParser.MODEL);
                writeNodeHeightPriorModelRef(prior, writer);
                writer.writeCloseTag(CoalescentLikelihoodParser.MODEL);
                writer.writeOpenTag(CoalescentLikelihoodParser.POPULATION_TREE);
                writer.writeIDref(DefaultTreeModel.TREE_MODEL, prefix + DefaultTreeModel.TREE_MODEL);
                writer.writeCloseTag(CoalescentLikelihoodParser.POPULATION_TREE);
                writer.writeCloseTag(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD);
        }
    }

    void writeNodeHeightPriorModelRef(PartitionTreePrior prior, XMLWriter writer) {
        TreePriorType treePrior = prior.getNodeHeightPrior();
        String priorPrefix = prior.getPrefix();

        switch (treePrior) {
            case CONSTANT:
                writer.writeIDref(ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL, priorPrefix + "constant");
                break;
            case EXPONENTIAL:
                writer.writeIDref(ExponentialGrowthModelParser.EXPONENTIAL_GROWTH_MODEL, priorPrefix + "exponential");
                break;
            case LOGISTIC:
                writer.writeIDref(LogisticGrowthModelParser.LOGISTIC_GROWTH_MODEL, priorPrefix + "logistic");
                break;
            case EXPANSION:
                writer.writeIDref(ExpansionModelParser.EXPANSION_MODEL, priorPrefix + "expansion");
                break;
//            case SKYLINE:
//                writer.writeIDref(BayesianSkylineLikelihoodParser.SKYLINE_LIKELIHOOD, priorPrefix + "skyline");
//                break;
            case GMRF_SKYRIDE:
                writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD, priorPrefix + "skyride");
                break;
            case SKYGRID:
            case SKYGRID_HMC:
                writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYGRID_LIKELIHOOD, priorPrefix + "skygrid");
                break;
            case YULE:
            case YULE_CALIBRATION:
                writer.writeIDref(YuleModelParser.YULE_MODEL, priorPrefix + YuleModelParser.YULE);
                break;
            case BIRTH_DEATH:
            case BIRTH_DEATH_INCOMPLETE_SAMPLING:
                writer.writeIDref(BirthDeathModelParser.BIRTH_DEATH_MODEL, priorPrefix + BirthDeathModelParser.BIRTH_DEATH);
                break;
            case BIRTH_DEATH_SERIAL_SAMPLING:
                writer.writeIDref(BirthDeathSerialSamplingModelParser.BIRTH_DEATH_SERIAL_MODEL,
                        priorPrefix + BirthDeathSerialSamplingModelParser.BDSS);
                break;
//            case BIRTH_DEATH_BASIC_REPRODUCTIVE_NUMBER:
//                writer.writeIDref(BirthDeathEpidemiologyModelParser.BIRTH_DEATH_EPIDEMIOLOGY,
//                        priorPrefix + BirthDeathEpidemiologyModelParser.BIRTH_DEATH_EPIDEMIOLOGY);
//                break;
            default:
                throw new IllegalArgumentException("No tree prior has been specified so cannot refer to it");
        }
    }

    void writeMultiLociTreePriors(PartitionTreePrior prior, XMLWriter writer) {

        String priorPrefix = prior.getPrefix();

        if (prior.getNodeHeightPrior() == TreePriorType.SKYGRID ||
                prior.getNodeHeightPrior() == TreePriorType.SKYGRID_HMC) {

            writer.writeComment("Generate a gmrfSkyGridLikelihood for the Bayesian SkyGrid process");
            writer.writeOpenTag(
                    GMRFSkyrideLikelihoodParser.SKYGRID_LIKELIHOOD,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, priorPrefix + "skygrid"),
                    }
            );

            int skyGridIntervalCount = prior.getSkyGridCount();
            double skyGridInterval = prior.getSkyGridInterval();

            writer.writeOpenTag(GMRFSkyrideLikelihoodParser.POPULATION_PARAMETER);
            writer.writeComment("skygrid.logPopSize is in log units unlike other popSize");
            writeParameter(prior.getParameter(GMRFSkyrideLikelihoodParser.SKYGRID_LOGPOPSIZE), skyGridIntervalCount, writer);
            writer.writeCloseTag(GMRFSkyrideLikelihoodParser.POPULATION_PARAMETER);

            writer.writeOpenTag(GMRFSkyrideLikelihoodParser.PRECISION_PARAMETER);
            writeParameter(prior.getParameter(GMRFSkyrideLikelihoodParser.SKYGRID_PRECISION), 1, writer);
            writer.writeCloseTag(GMRFSkyrideLikelihoodParser.PRECISION_PARAMETER);

            writer.writeOpenTag(GMRFSkyrideLikelihoodParser.NUM_GRID_POINTS);
            Parameter numGridPoint = prior.getParameter("skygrid.numGridPoints");
            numGridPoint.setInitial(skyGridIntervalCount - 1);
            writeParameter(numGridPoint, 1, writer);
            writer.writeCloseTag(GMRFSkyrideLikelihoodParser.NUM_GRID_POINTS);

            writer.writeOpenTag(GMRFSkyrideLikelihoodParser.CUT_OFF);
            Parameter cutOff = prior.getParameter("skygrid.cutOff");
            cutOff.setInitial(skyGridInterval);
            writeParameter(cutOff, 1, writer);
            writer.writeCloseTag(GMRFSkyrideLikelihoodParser.CUT_OFF);

            writer.writeOpenTag(GMRFSkyrideLikelihoodParser.POPULATION_TREE);
            // TODO Add all linked trees
            if (options.isShareSameTreePrior()) {
                for (PartitionTreeModel thisModel : options.getPartitionTreeModels()) {
                    writer.writeIDref(DefaultTreeModel.TREE_MODEL, thisModel.getPrefix() + DefaultTreeModel.TREE_MODEL);
                }
            } else {
                writer.writeIDref(DefaultTreeModel.TREE_MODEL, options.getPartitionTreeModels(prior).get(0).getPrefix() + DefaultTreeModel.TREE_MODEL);
            }
            writer.writeCloseTag(GMRFSkyrideLikelihoodParser.POPULATION_TREE);

            writer.writeCloseTag(GMRFSkyrideLikelihoodParser.SKYGRID_LIKELIHOOD);

            //writing the gamma prior here so will need to prevent another one from being written in the priors block
            //key use: using HMC on the skygrid parameters

            Parameter parameter = prior.getParameter(GMRFSkyrideLikelihoodParser.SKYGRID_PRECISION);
            writer.writeOpenTag(PriorParsers.GAMMA_PRIOR,
                    new Attribute[]{
                            new Attribute.Default<>(XMLParser.ID, GMRFSkyrideLikelihoodParser.SKYGRID_PRECISION_PRIOR),
                            new Attribute.Default<>(GammaDistributionModelParser.SHAPE, parameter.shape),
                            new Attribute.Default<>(GammaDistributionModelParser.SCALE, parameter.scale),
                            new Attribute.Default<>(GammaDistributionModelParser.OFFSET, parameter.offset)
                    }
            );
            writer.writeIDref(ParameterParser.PARAMETER, GMRFSkyrideLikelihoodParser.SKYGRID_PRECISION);
            writer.writeCloseTag(PriorParsers.GAMMA_PRIOR);

            //add gradient information to XML file in case of an HMC transition kernel mix
            if (prior.getNodeHeightPrior() == TreePriorType.SKYGRID_HMC) {

                writer.writeOpenTag(GMRFSkyrideGradientParser.NAME,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, "gmrfGradientPop"),
                                new Attribute.Default<String>(GMRFSkyrideGradientParser.WRT_PARAMETER, "logPopulationSizes")
                        }
                );
                writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYGRID_LIKELIHOOD, "skygrid");
                writer.writeCloseTag(GMRFSkyrideGradientParser.NAME);

                writer.writeOpenTag(CompoundParameterParser.COMPOUND_PARAMETER,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, "skygrid.parameters")
                        }
                );
                writer.writeIDref(ParameterParser.PARAMETER, GMRFSkyrideLikelihoodParser.SKYGRID_PRECISION);
                writer.writeIDref(ParameterParser.PARAMETER, GMRFSkyrideLikelihoodParser.SKYGRID_LOGPOPSIZE);
                writer.writeCloseTag(CompoundParameterParser.COMPOUND_PARAMETER);

                writer.writeOpenTag(GMRFSkyrideGradientParser.NAME,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, "gmrfGradientPrec"),
                                new Attribute.Default<String>(GMRFSkyrideGradientParser.WRT_PARAMETER, "precision")
                        }
                );
                writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYGRID_LIKELIHOOD, "skygrid");
                writer.writeCloseTag(GMRFSkyrideGradientParser.NAME);

                writer.writeOpenTag(JointGradientParser.JOINT_GRADIENT,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, "joint.skygrid.precision")
                        }
                );
                writer.writeIDref(GMRFSkyrideGradientParser.NAME, "gmrfGradientPrec");
                writer.writeOpenTag(GradientWrapperParser.NAME);
                writer.writeIDref(PriorParsers.GAMMA_PRIOR, GMRFSkyrideLikelihoodParser.SKYGRID_PRECISION_PRIOR);
                writer.writeIDref(ParameterParser.PARAMETER, GMRFSkyrideLikelihoodParser.SKYGRID_PRECISION);
                writer.writeCloseTag(GradientWrapperParser.NAME);
                writer.writeCloseTag(JointGradientParser.JOINT_GRADIENT);

                writer.writeOpenTag(CompoundGradientParser.COMPOUND_GRADIENT,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID,  "full.skygrid.gradient")
                        }
                );
                writer.writeIDref(JointGradientParser.JOINT_GRADIENT, "joint.skygrid.precision");
                writer.writeIDref(GMRFSkyrideGradientParser.NAME, "gmrfGradientPop");
                writer.writeCloseTag(CompoundGradientParser.COMPOUND_GRADIENT);

            }
        }
    }

    void writeParameterLog(PartitionTreePrior prior, XMLWriter writer) {

        String priorPrefix = prior.getPrefix();

        switch (prior.getNodeHeightPrior()) {

            case CONSTANT:
                writeParameterRef(priorPrefix + "constant.popSize", writer);
                break;
            case EXPONENTIAL:
                writeParameterRef(priorPrefix + "exponential.popSize", writer);
                if (prior.getParameterization() == TreePriorParameterizationType.GROWTH_RATE) {
                    writeParameterRef(priorPrefix + "exponential.growthRate", writer);
                } else {
                    writeParameterRef(priorPrefix + "exponential.doublingTime", writer);
                }
                break;
            case LOGISTIC:
                writeParameterRef(priorPrefix + "logistic.popSize", writer);
                if (prior.getParameterization() == TreePriorParameterizationType.GROWTH_RATE) {
                    writeParameterRef(priorPrefix + "logistic.growthRate", writer);
                } else {
                    writeParameterRef(priorPrefix + "logistic.doublingTime", writer);
                }
                writeParameterRef(priorPrefix + "logistic.t50", writer);
                break;
            case EXPANSION:
                writeParameterRef(priorPrefix + "expansion.popSize", writer);
                if (prior.getParameterization() == TreePriorParameterizationType.GROWTH_RATE) {
                    writeParameterRef(priorPrefix + "expansion.growthRate", writer);
                } else {
                    writeParameterRef(priorPrefix + "expansion.doublingTime", writer);
                }
                writeParameterRef(priorPrefix + "expansion.ancestralProportion", writer);
                break;
//            case SKYLINE:
//                writeParameterRef(priorPrefix + "skyline.popSize", writer);
//                writeParameterRef(priorPrefix + "skyline.groupSize", writer);
//                break;
            case SKYGRID:
            case SKYGRID_HMC:
                writeParameterRef(priorPrefix + GMRFSkyrideLikelihoodParser.SKYGRID_PRECISION, writer);
                writeParameterRef(priorPrefix + GMRFSkyrideLikelihoodParser.SKYGRID_LOGPOPSIZE, writer);
                writeParameterRef(priorPrefix + "skygrid.cutOff", writer);
//                writeParameterRef(priorPrefix + "skygrid.groupSize", writer);
                break;
            case GMRF_SKYRIDE:
                writeParameterRef(priorPrefix + "skyride.precision", writer);
                writeParameterRef(priorPrefix + "skyride.logPopSize", writer);
                writeParameterRef(priorPrefix + "skyride.groupSize", writer);
                break;
            case YULE:
            case YULE_CALIBRATION:
                writeParameterRef(priorPrefix + "yule.birthRate", writer);
                break;
            case BIRTH_DEATH:
            case BIRTH_DEATH_INCOMPLETE_SAMPLING:
                writeParameterRef(priorPrefix + BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME, writer);
                writeParameterRef(priorPrefix + BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME, writer);
                if (prior.getNodeHeightPrior() == TreePriorType.BIRTH_DEATH_INCOMPLETE_SAMPLING)
                    writeParameterRef(priorPrefix + BirthDeathModelParser.BIRTH_DEATH + "."
                            + BirthDeathModelParser.SAMPLE_PROB, writer);
                break;
            case BIRTH_DEATH_SERIAL_SAMPLING:
                writeParameterRef(priorPrefix + BirthDeathSerialSamplingModelParser.BDSS + "."
                        + BirthDeathSerialSamplingModelParser.LAMBDA, writer);
                writeParameterRef(priorPrefix + BirthDeathSerialSamplingModelParser.BDSS + "."
                        + BirthDeathSerialSamplingModelParser.RELATIVE_MU, writer);
                //Issue 656
//                writeParameterRef(modelPrefix + BirthDeathSerialSamplingModelParser.BDSS + "."
//                        + BirthDeathSerialSamplingModelParser.SAMPLE_PROBABILITY, writer);
                writeParameterRef(priorPrefix + BirthDeathSerialSamplingModelParser.BDSS + "."
                        + BirthDeathSerialSamplingModelParser.PSI, writer);
                writeParameterRef(priorPrefix + BirthDeathSerialSamplingModelParser.BDSS + "."
                        + BirthDeathSerialSamplingModelParser.ORIGIN, writer);
                break;
//            case BIRTH_DEATH_BASIC_REPRODUCTIVE_NUMBER:
//                writeParameterRef(priorPrefix + BirthDeathEpidemiologyModelParser.R0, writer);
//                writeParameterRef(priorPrefix + BirthDeathEpidemiologyModelParser.RECOVERY_RATE, writer);
//                writeParameterRef(priorPrefix + BirthDeathEpidemiologyModelParser.SAMPLING_PROBABILITY, writer);
//                writeParameterRef(priorPrefix + BirthDeathEpidemiologyModelParser.ORIGIN, writer);
            default:
                throw new IllegalArgumentException("No tree prior has been specified so cannot refer to it");
        }

    }

    public static void writePriorLikelihoodReferenceLog(PartitionTreePrior prior, PartitionTreeModel model, XMLWriter writer) {

        String prefix = model.getPrefix();

        switch (prior.getNodeHeightPrior()) {

            case YULE:
            case YULE_CALIBRATION:
            case BIRTH_DEATH:
            case BIRTH_DEATH_INCOMPLETE_SAMPLING:
            case BIRTH_DEATH_SERIAL_SAMPLING:
//            case BIRTH_DEATH_BASIC_REPRODUCTIVE_NUMBER:
                writer.writeIDref(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD, prefix + "speciation");
                break;
//            case SKYLINE:
//                writer.writeIDref(BayesianSkylineLikelihoodParser.SKYLINE_LIKELIHOOD, prefix + "skyline");
////                writer.writeIDref(ExponentialMarkovModel.EXPONENTIAL_MARKOV_MODEL, modelPrefix + "eml1");
//                break;
            case GMRF_SKYRIDE:
                writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD, prefix + "skyride");
                break;
            case SKYGRID:
            case SKYGRID_HMC:
                writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD, prefix + "skygrid");
                // only 1 coalescent, so write it separately after this method
                break;
            case LOGISTIC:
//                writer.writeIDref(BooleanLikelihoodParser.BOOLEAN_LIKELIHOOD, modelPrefix + "booleanLikelihood1");
                writer.writeIDref(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, prefix + COALESCENT);
                break;
            default:
                writer.writeIDref(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, prefix + COALESCENT);
        }
    }

    // id is written in writePriorLikelihood (PartitionTreePrior prior, PartitionTreeModel model, XMLWriter writer)
    public void writePriorLikelihoodReference(PartitionTreePrior prior, PartitionTreeModel model, XMLWriter writer) {

        String prefix = model.getPrefix();

        switch (prior.getNodeHeightPrior()) {

            case YULE:
            case YULE_CALIBRATION:
            case BIRTH_DEATH:
            case BIRTH_DEATH_INCOMPLETE_SAMPLING:
            case BIRTH_DEATH_SERIAL_SAMPLING:
//            case BIRTH_DEATH_BASIC_REPRODUCTIVE_NUMBER:
                writer.writeIDref(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD, prefix + "speciation");
                break;
//            case SKYLINE:
//                writer.writeIDref(BayesianSkylineLikelihoodParser.SKYLINE_LIKELIHOOD, prefix + "skyline");
//                writer.writeIDref(ExponentialMarkovModel.EXPONENTIAL_MARKOV_MODEL, prefix + "eml1");
//                break;
            case GMRF_SKYRIDE:
                writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD, prefix + "skyride");
                break;
            case SKYGRID:
            case SKYGRID_HMC:
//                writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD, prefix + "skygrid");
                // only 1 coalescent, so write it separately after this method
                break;
            default:
                writer.writeIDref(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, prefix + COALESCENT);
        }
    }

    public void writeMultiLociLikelihoodReference(PartitionTreePrior prior, XMLWriter writer) {
        String prefix = prior.getPrefix();
        if ((prior.getNodeHeightPrior() == TreePriorType.SKYGRID || prior.getNodeHeightPrior() == TreePriorType.SKYGRID_HMC)) {
            writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYGRID_LIKELIHOOD, prefix + "skygrid");
        }

    }
}
