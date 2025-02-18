/*
 * PartitionSubstitutionModel.java
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

import dr.app.beauti.components.continuous.ContinuousModelExtensionType;
import dr.evomodel.substmodel.aminoacid.AminoAcidModelType;
import dr.evomodel.substmodel.nucleotide.NucModelType;
import dr.app.beauti.components.continuous.ContinuousSubstModelType;
import dr.app.beauti.components.discrete.DiscreteSubstModelType;
import dr.app.beauti.types.*;
import dr.evolution.datatype.AminoAcids;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.stats.Variate;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class PartitionSubstitutionModel extends PartitionOptions {
    private static final long serialVersionUID = -2570346396317131108L;

    // Instance variables

    public static final String GTR_RATES = "gtr.rates";
    public static final String[] GTR_RATE_NAMES = {"gtr.AC", "gtr.AG", "gtr.AT", "gtr.CG", "gtr.GT"};
    private static final String[] GTR_TRANSITIONS = {"A-C", "A-G", "A-T", "C-G", "G-T"};

    private NucModelType nucSubstitutionModel = NucModelType.HKY;
    private AminoAcidModelType aaSubstitutionModel = AminoAcidModelType.BLOSUM_62;
    private BinaryModelType binarySubstitutionModel = BinaryModelType.BIN_SIMPLE;
    private DiscreteSubstModelType discreteSubstType = DiscreteSubstModelType.SYM_SUBST;
    private ContinuousSubstModelType continuousSubstModelType = ContinuousSubstModelType.HOMOGENOUS;
    private ContinuousModelExtensionType continuousExtensionType = ContinuousModelExtensionType.NONE;

    private final int dataType;
    private int continuousTraitCount;
    private final int extendedTraitCount;

    private final TraitData traitData;

    private boolean activateBSSVS = false;
    public boolean useAmbiguitiesTreeLikelihood = false;

    private FrequencyPolicyType frequencyPolicy = FrequencyPolicyType.ESTIMATED;
    private boolean gammaHetero = false;
    private int gammaCategories = 4;
    private boolean invarHetero = false;
    private boolean equalWeights = false;
    private String codonHeteroPattern = null;
    private boolean unlinkedSubstitutionModel = true;
    private boolean unlinkedHeterogeneityModel = true;
    private boolean unlinkedFrequencyModel = true;

    private boolean dolloModel = false;

    private boolean isLatitudeLongitude = false;
    private boolean isIndependent = false;
    private double jitterWindow = 0.0;

    public TraitData getTraitData() {
        return traitData;
    }

    public PartitionSubstitutionModel(BeautiOptions options, String name, AbstractPartitionData partition) {
        super(options, name);

        if (partition.getTraits() != null && partition.getDataType().getType() == DataType.CONTINUOUS) {
            continuousTraitCount = partition.getTraits().size();
        } else {
            continuousTraitCount = 0;
        }

        extendedTraitCount = continuousTraitCount;

        if (partition.getTraits() != null && partition.getDataType().getType() == DataType.GENERAL) {
            traitData = partition.getTraits().get(0);
        } else {
            traitData = null;
        }

        dataType = partition.getDataType().getType();
        initModelParametersAndOpererators();
    }

    /**
     * A copy constructor
     *
     * @param options the beauti options
     * @param name    the name of the new model
     * @param source  the source model
     */
    public PartitionSubstitutionModel(BeautiOptions options, String name, PartitionSubstitutionModel source) {
        super(options, name);

        nucSubstitutionModel = source.nucSubstitutionModel;
        aaSubstitutionModel = source.aaSubstitutionModel;
        binarySubstitutionModel = source.binarySubstitutionModel;
        discreteSubstType = source.discreteSubstType;
        continuousSubstModelType = source.continuousSubstModelType;
        continuousExtensionType = source.continuousExtensionType;

        continuousTraitCount = source.continuousTraitCount;
        extendedTraitCount = source.extendedTraitCount;
        dataType = source.dataType;

        traitData = source.traitData;

        activateBSSVS = source.activateBSSVS;
        useAmbiguitiesTreeLikelihood = source.useAmbiguitiesTreeLikelihood;

        frequencyPolicy = source.frequencyPolicy;
        gammaHetero = source.gammaHetero;
        gammaCategories = source.gammaCategories;
        equalWeights = source.equalWeights;
        invarHetero = source.invarHetero;
        codonHeteroPattern = source.codonHeteroPattern;
        unlinkedSubstitutionModel = source.unlinkedSubstitutionModel;
        unlinkedHeterogeneityModel = source.unlinkedHeterogeneityModel;
        unlinkedFrequencyModel = source.unlinkedFrequencyModel;

        dolloModel = source.dolloModel;

        initModelParametersAndOpererators();
    }

    public PartitionSubstitutionModel(BeautiOptions options, String name) {
        super(options, name);
        continuousTraitCount = 0;
        extendedTraitCount = continuousTraitCount;

        traitData = null;

        dataType = DataType.NUCLEOTIDES;

        initModelParametersAndOpererators();
    }

    // only init in PartitionSubstitutionModel
    @Override
    public void initModelParametersAndOpererators() {
        double substWeights = 1.0;
        if (dataType == DataType.NUCLEOTIDES || dataType == DataType.AMINO_ACIDS) {
            //Substitution model parameters
            if (options.useNewFrequenciesPrior()) {
                createNonNegativeParameterDirichletPrior("frequencies", "base frequencies", this, 4, 1.0, true);
                createNonNegativeParameterDirichletPrior("CP1.frequencies", "base frequencies for codon position 1", this, 4, 1.0, true);
                createNonNegativeParameterDirichletPrior("CP2.frequencies", "base frequencies for codon position 2", this, 4, 1.0, true);
                createNonNegativeParameterDirichletPrior("CP1+2.frequencies", "base frequencies for codon positions 1 & 2", this, 4, 1.0, true);
                createNonNegativeParameterDirichletPrior("CP3.frequencies", "base frequencies for codon position 3", this, 4, 1.0, true);
            } else {
                createZeroOneParameterUniformPrior("frequencies", "base frequencies", 0.25, true);
                createZeroOneParameterUniformPrior("CP1.frequencies", "base frequencies for codon position 1", 0.25, true);
                createZeroOneParameterUniformPrior("CP2.frequencies", "base frequencies for codon position 2", 0.25, true);
                createZeroOneParameterUniformPrior("CP1+2.frequencies", "base frequencies for codon positions 1 & 2", 0.25, true);
                createZeroOneParameterUniformPrior("CP3.frequencies", "base frequencies for codon position 3", 0.25, true);
            }

            //This prior is moderately diffuse with a median of 2.718
            createParameterLognormalPrior("kappa", "HKY transition-transversion parameter",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 2.0, 1.0, 1.25, 0.0, true);
            createParameterLognormalPrior("CP1.kappa", "HKY transition-transversion parameter for codon position 1",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 2.0, 1.0, 1.25, 0.0, true);
            createParameterLognormalPrior("CP2.kappa", "HKY transition-transversion parameter for codon position 2",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 2.0, 1.0, 1.25, 0.0, true);
            createParameterLognormalPrior("CP1+2.kappa", "HKY transition-transversion parameter for codon positions 1 & 2",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 2.0, 1.0, 1.25, 0.0, true);
            createParameterLognormalPrior("CP3.kappa", "HKY transition-transversion parameter for codon position 3",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 2.0, 1.0, 1.25, 0.0, true);

            createParameterLognormalPrior("kappa1", "TN93 1st transition-transversion parameter",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 2.0, 1.0, 1.25, 0.0, true);
            createParameterLognormalPrior("CP1.kappa1", "TN93 1st transition-transversion parameter for codon position 1",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 2.0, 1.0, 1.25, 0.0, true);
            createParameterLognormalPrior("CP2.kappa1", "TN93 1st transition-transversion parameter for codon position 2",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 2.0, 1.0, 1.25, 0.0, true);
            createParameterLognormalPrior("CP1+2.kappa1", "TN93 1st transition-transversion parameter for codon positions 1 & 2",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 2.0, 1.0, 1.25, 0.0, true);
            createParameterLognormalPrior("CP3.kappa1", "TN93 1st transition-transversion parameter for codon position 3",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 2.0, 1.0, 1.25, 0.0, true);

            createParameterLognormalPrior("kappa2", "TN93 2nd transition-transversion parameter",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 2.0, 1.0, 1.25, 0.0, true);
            createParameterLognormalPrior("CP1.kappa2", "TN93 2nd transition-transversion parameter for codon position 1",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 2.0, 1.0, 1.25, 0.0, true);
            createParameterLognormalPrior("CP2.kappa2", "TN93 2nd transition-transversion parameter for codon position 2",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 2.0, 1.0, 1.25, 0.0, true);
            createParameterLognormalPrior("CP1+2.kappa2", "TN93 2nd transition-transversion parameter for codon positions 1 & 2",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 2.0, 1.0, 1.25, 0.0, true);
            createParameterLognormalPrior("CP3.kappa2", "TN93 2nd transition-transversion parameter for codon position 3",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 2.0, 1.0, 1.25, 0.0, true);

//        createParameter("frequencies", "GTR base frequencies", UNITY_SCALE, 0.25, 0.0, 1.0);
//        createParameter("CP1.frequencies", "GTR base frequencies for codon position 1", UNITY_SCALE, 0.25, 0.0, 1.0);
//        createParameter("CP2.frequencies", "GTR base frequencies for codon position 2", UNITY_SCALE, 0.25, 0.0, 1.0);
//        createParameter("CP1+2.frequencies", "GTR base frequencies for codon positions 1 & 2", UNITY_SCALE, 0.25, 0.0, 1.0);
//        createParameter("CP3.frequencies", "GTR base frequencies for codon position 3", UNITY_SCALE, 0.25, 0.0, 1.0);

            // create the relative rate parameters for the GTR rate matrix

            // if (!options.classicOperatorsAndPriors && options.NEW_GTR_PARAMETERIZATION) {
            createNonNegativeParameterDirichletPrior(GTR_RATES, "GTR transition rates parameter",
                    this, 6, 6.0, true);
            for (int i = 1; i <= 3; i++) {
                createNonNegativeParameterDirichletPrior("CP" + i + "." + GTR_RATES, "GTR transition rates parameter",
                        this, 6, 6.0, true);
            }
            createNonNegativeParameterDirichletPrior("CP1+2." + GTR_RATES, "GTR transition rates parameter",
                    this, 6, 6.0, true);

            //} else {
            for (int j = 0; j < 5; j++) {
                if (j == 1) { // ag
                    createParameterGammaPrior(GTR_RATE_NAMES[j], "GTR " + GTR_TRANSITIONS[j] + " substitution parameter",
                            PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.05, 20, false, true);
                } else {
                    createParameterGammaPrior(GTR_RATE_NAMES[j], "GTR " + GTR_TRANSITIONS[j] + " substitution parameter",
                            PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.05, 10, false, true);
                }

                for (int i = 1; i <= 3; i++) {
                    if (j == 1) { // ag
                        createParameterGammaPrior("CP" + i + "." + GTR_RATE_NAMES[j], "GTR " + GTR_TRANSITIONS[j] + " substitution parameter for codon position " + i,
                                PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.05, 20, false, true);
                    } else {
                        createParameterGammaPrior("CP" + i + "." + GTR_RATE_NAMES[j], "GTR " + GTR_TRANSITIONS[j] + " substitution parameter for codon position " + i,
                                PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.05, 10, false, true);
                    }
                }

                if (j == 1) { // ag
                    createParameterGammaPrior("CP1+2." + GTR_RATE_NAMES[j], "GTR " + GTR_TRANSITIONS[j] + " substitution parameter for codon positions 1 & 2",
                            PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.05, 20, false, true);
                } else {
                    createParameterGammaPrior("CP1+2." + GTR_RATE_NAMES[j], "GTR " + GTR_TRANSITIONS[j] + " substitution parameter for codon positions 1 & 2",
                            PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.05, 10, false, true);
                }
            }
            //}

//        createParameter("frequencies", "Binary Simple frequencies", UNITY_SCALE, 0.5, 0.0, 1.0);
//
//        createParameter("frequencies", "Binary Covarion frequencies of the visible states", UNITY_SCALE, 0.5, 0.0, 1.0);
            createZeroOneParameterUniformPrior("hfrequencies", "Binary Covarion frequencies of the hidden rates", 0.5, true);
            createZeroOneParameterUniformPrior("bcov.alpha", "Binary Covarion rate of evolution in slow mode", 0.5, true);
            createParameterGammaPrior("bcov.s", "Binary Covarion rate of flipping between slow and fast modes",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.05, 10, false, true);

            createParameterExponentialPrior("alpha", "gamma shape parameter",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.5, 0.0, true);
            createParameterExponentialPrior("CP1.alpha", "gamma shape parameter for codon position 1",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.5, 0.0, true);
            createParameterExponentialPrior("CP2.alpha", "gamma shape parameter for codon position 2",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.5, 0.0, true);
            createParameterExponentialPrior("CP1+2.alpha", "gamma shape parameter for codon positions 1 & 2",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.5, 0.0, true);
            createParameterExponentialPrior("CP3.alpha", "gamma shape parameter for codon position 3",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.5, 0.0, true);

            createZeroOneParameterUniformPrior("pInv", "proportion of invariant sites parameter", 0.5, true);
            createZeroOneParameterUniformPrior("CP1.pInv", "proportion of invariant sites parameter for codon position 1", 0.5, true);
            createZeroOneParameterUniformPrior("CP2.pInv", "proportion of invariant sites parameter for codon position 2", 0.5, true);
            createZeroOneParameterUniformPrior("CP1+2.pInv", "proportion of invariant sites parameter for codon positions 1 & 2", 0.5, true);
            createZeroOneParameterUniformPrior("CP3.pInv", "proportion of invariant sites parameter for codon position 3", 0.5, true);

            // nu parameters are an alternative parameterization of relative rates suitable for using with a Dirichlet prior.
            // They sum to 1 and their product with the evolutionary rate is weighted by the partitions size.
            createZeroOneParameterUniformPrior("nu", "relative rate parameter", 1.0, true);
            createZeroOneParameterUniformPrior("CP1.nu", "relative rate parameter for codon position 1", 1.0, true);
            createZeroOneParameterUniformPrior("CP2.nu", "relative rate parameter for codon position 2", 1.0, true);
            createZeroOneParameterUniformPrior("CP1+2.nu", "relative rate parameter for codon positions 1 & 2", 1.0, true);
            createZeroOneParameterUniformPrior("CP3.nu", "relative rate parameter for codon position 3", 1.0, true);

            createNonNegativeParameterInfinitePrior("mu", "relative rate parameter", PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, true);
            createNonNegativeParameterInfinitePrior("CP1.mu", "relative rate parameter for codon position 1",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, true);
            createNonNegativeParameterInfinitePrior("CP2.mu", "relative rate parameter for codon position 2",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, true);
            createNonNegativeParameterInfinitePrior("CP1+2.mu", "relative rate parameter for codon positions 1 & 2",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, true);
            createNonNegativeParameterInfinitePrior("CP3.mu", "relative rate parameter for codon position 3",
                    PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, true);

            createScaleOperator("kappa", demoTuning, substWeights);
            createScaleOperator("CP1.kappa", demoTuning, substWeights);
            createScaleOperator("CP2.kappa", demoTuning, substWeights);
            createScaleOperator("CP1+2.kappa", demoTuning, substWeights);
            createScaleOperator("CP3.kappa", demoTuning, substWeights);

            createScaleOperator("kappa1", demoTuning, substWeights);
            createScaleOperator("CP1.kappa1", demoTuning, substWeights);
            createScaleOperator("CP2.kappa1", demoTuning, substWeights);
            createScaleOperator("CP1+2.kappa1", demoTuning, substWeights);
            createScaleOperator("CP3.kappa1", demoTuning, substWeights);

            createScaleOperator("kappa2", demoTuning, substWeights);
            createScaleOperator("CP1.kappa2", demoTuning, substWeights);
            createScaleOperator("CP2.kappa2", demoTuning, substWeights);
            createScaleOperator("CP1+2.kappa2", demoTuning, substWeights);
            createScaleOperator("CP3.kappa2", demoTuning, substWeights);

            createOperator("frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
            createOperator("CP1.frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
            createOperator("CP2.frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
            createOperator("CP1+2.frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
            createOperator("CP3.frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);

            // if (!options.classicOperatorsAndPriors && options.NEW_GTR_PARAMETERIZATION) {
            createOperator("deltaGTR", "gtr.rates",
                    "Change GTR transition rates relative to each other maintaining mean",
                    GTR_RATES,
                    OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
            for (int j = 1; j <= 3; j++) {
                createOperator("CP" + j + ".deltaGTR", "CP" + j + ".gtr.rates",
                        "Change GTR transition rates relative to each other maintaining mean",
                        "CP" + j + "." + GTR_RATES,
                        OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
            }
            createOperator("CP1+2.deltaGTR", "CP1+2.gtr.rates",
                    "Change GTR transition rates relative to each other maintaining mean",
                    "CP1+2." + GTR_RATES,
                    OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
            // } else {
            for (String rateName : GTR_RATE_NAMES) {
                createScaleOperator(rateName, demoTuning, substWeights);
                for (int j = 1; j <= 3; j++) {
                    createScaleOperator("CP" + j + "." + rateName, demoTuning, substWeights);
                }
                createScaleOperator("CP1+2." + rateName, demoTuning, substWeights);
            }
            // }

            createScaleOperator("alpha", demoTuning, substWeights);
            for (int i = 1; i <= 3; i++) {
                createScaleOperator("CP" + i + ".alpha", demoTuning, substWeights);
            }
            createScaleOperator("CP1+2.alpha", demoTuning, substWeights);

            // if (!options.classicOperatorsAndPriors && LOGIT_PINV_KERNEL) { // a switch at the top of BeautiOptions
            createOperator("rwPInv", "rwPInv", "Random walk on pInv in logit space", "pInv", OperatorType.RANDOM_WALK_LOGIT, demoTuning, substWeights);
            for (int i = 1; i <= 3; i++) {
                createOperator("CP" + i + ".rwPInv", "CP" + i + ".rwPInv", "Random walk on pInv in logit space", "CP" + i + ".pInv", OperatorType.RANDOM_WALK_LOGIT, demoTuning, substWeights);
            }
            createOperator("CP1+2.rwPInv", "CP1+2.rwPInv", "Random walk on pInv in logit space", "CP1+2.pInv", OperatorType.RANDOM_WALK_LOGIT, demoTuning, substWeights);
            // } else {
            // old (and not very appropriate scale operator)
            createOperator("CP1+2.uniformPInv", "CP1+2.uniformPInv", "Random walk on pInv in logit space", "pInv", OperatorType.UNIFORM, demoTuning, substWeights);
            for (int i = 1; i <= 3; i++) {
                createOperator("CP" + i + ".uniformPInv", "CP" + i + ".uniformPInv", "Random walk on pInv in logit space", "CP" + i + ".pInv", OperatorType.UNIFORM, demoTuning, substWeights);
            }
            createOperator("CP1+2.uniformPInv", "CP1+2.uniformPInv", "Random walk on pInv in logit space", "CP1+2.pInv", OperatorType.UNIFORM, demoTuning, substWeights);
            // }
        } else if (dataType == DataType.COVARION) {
            createZeroOneParameterUniformPrior("nu", "relative rate parameter", 1.0, true);
            createNonNegativeParameterInfinitePrior("mu", "relative rate parameter", PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, true);
            createScaleOperator("bcov.alpha", demoTuning, substWeights);
            createScaleOperator("bcov.s", demoTuning, substWeights);
            createOperator("hfrequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
        } else {
            createZeroOneParameterUniformPrior("nu", "relative rate parameter", 1.0, true);
            createNonNegativeParameterInfinitePrior("mu", "relative rate parameter", PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, true);
        }
    }

    ////////////////////////////////////////////////////////////////

    @Override
    public List<Parameter> selectParameters(List<Parameter> params) {
//        setAvgRootAndRate();
        boolean includeRelativeRates = getCodonPartitionCount() > 1;//TODO check
        switch (getDataType().getType()) {
            case DataType.NUCLEOTIDES:
                if (includeRelativeRates && unlinkedSubstitutionModel) {
                    if (codonHeteroPattern.equals("123")) {
                        for (int i = 1; i <= getCodonPartitionCount(); i++) {
                            switch (nucSubstitutionModel) {
                                case JC:
                                    break;
                                case HKY:
                                    params.add(getParameter("CP" + i + ".kappa"));
                                    break;
                                case TN93:
                                    params.add(getParameter("CP" + i + ".kappa1"));
                                    params.add(getParameter("CP" + i + ".kappa2"));
                                    break;
                                case GTR:
                                    if (options.useNewGTR()) {
                                        params.add(getParameter("CP" + i + "." + GTR_RATES));
                                    } else {
                                        for (String rateName : GTR_RATE_NAMES) {
                                            params.add(getParameter("CP" + i + "." + rateName));
                                        }
                                    }
                                    break;

                                default:
                                    throw new IllegalArgumentException("Unknown nucleotides substitution model");
                            }
                        }

                    } else if (codonHeteroPattern.equals("112")) {
                        switch (nucSubstitutionModel) {
                            case JC:
                                break;
                            case HKY:
                                params.add(getParameter("CP1+2.kappa"));
                                params.add(getParameter("CP3.kappa"));
                                break;
                            case TN93:
                                params.add(getParameter("CP1+2.kappa1"));
                                params.add(getParameter("CP3.kappa1"));
                                params.add(getParameter("CP1+2.kappa2"));
                                params.add(getParameter("CP3.kappa2"));
                                break;
                            case GTR:
                                if (options.useNewGTR()) {
                                    params.add(getParameter("CP1+2." + GTR_RATES));
                                    params.add(getParameter("CP3." + GTR_RATES));
                                } else {
                                    for (String rateName : GTR_RATE_NAMES) {
                                        params.add(getParameter("CP1+2." + rateName));
                                    }
                                    for (String rateName : GTR_RATE_NAMES) {
                                        params.add(getParameter("CP3." + rateName));
                                    }
                                }
                                break;

                            default:
                                throw new IllegalArgumentException("Unknown nucleotides substitution model");
                        }

                    } else {
                        throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                    }
                } else { // no codon partitioning, or unlinkedSubstitutionModel
                    switch (nucSubstitutionModel) {
                        case JC:
                            break;
                        case HKY:
                            params.add(getParameter("kappa"));
                            break;
                        case TN93:
                            params.add(getParameter("kappa1"));
                            params.add(getParameter("kappa2"));
                            break;
                        case GTR:
                            if (options.useNewGTR()) {
                                params.add(getParameter(GTR_RATES));
                            } else {

                                for (String rateName : GTR_RATE_NAMES) {
                                    params.add(getParameter(rateName));
                                }
                            }
                            break;

                        default:
                            throw new IllegalArgumentException("Unknown nucleotides substitution model");
                    }
                }

                // only AMINO_ACIDS not addFrequency
                addFrequencyParams(params, includeRelativeRates);
                break;

            case DataType.AMINO_ACIDS:
//                if (includeRelativeRates) {
//                    params.add(getParameter("mu"));
//                }
                break;

            case DataType.TWO_STATES:
            case DataType.COVARION:
                switch (binarySubstitutionModel) {
                    case BIN_SIMPLE:
                    case BIN_DOLLO:
                        break;

                    case BIN_COVARION:
//                        useAmbiguitiesTreeLikelihood = true;
                        params.add(getParameter("bcov.alpha"));
                        params.add(getParameter("bcov.s"));
                        params.add(getParameter("hfrequencies")); // no codon for binary
                        break;

                    default:
                        throw new IllegalArgumentException("Unknown binary substitution model");
                }
//                if (includeRelativeRates) {
//                    params.add(getParameter("mu"));
//                }

                // only AMINO_ACIDS not addFrequency
                addFrequencyParams(params, includeRelativeRates);
                break;

            case DataType.GENERAL:
                // This model is controlled by DiscreteTraitComponentOptions
                break;

            case DataType.CONTINUOUS:
                // This model is controlled by ContinuousTraitComponentOptions
                break;

            case DataType.DUMMY:
                //Do nothing
                break;

            case DataType.TREE:
                //Do nothing
                break;
            default:
                throw new IllegalArgumentException("Unknown data type");
        }

        // if gamma do shape move
        if (gammaHetero) {
            if (includeRelativeRates && unlinkedHeterogeneityModel) {
                if (codonHeteroPattern.equals("123")) {
                    params.add(getParameter("CP1.alpha"));
                    params.add(getParameter("CP2.alpha"));
                    params.add(getParameter("CP3.alpha"));
                } else if (codonHeteroPattern.equals("112")) {
                    params.add(getParameter("CP1+2.alpha"));
                    params.add(getParameter("CP3.alpha"));
                } else {
                    throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                }
            } else {
                params.add(getParameter("alpha"));
            }
        }
        // if pinv do pinv move
        if (invarHetero) {
            if (includeRelativeRates && unlinkedHeterogeneityModel) {
                if (codonHeteroPattern.equals("123")) {
                    params.add(getParameter("CP1.pInv"));
                    params.add(getParameter("CP2.pInv"));
                    params.add(getParameter("CP3.pInv"));
                } else if (codonHeteroPattern.equals("112")) {
                    params.add(getParameter("CP1+2.pInv"));
                    params.add(getParameter("CP3.pInv"));
                } else {
                    throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                }
            } else {
                params.add(getParameter("pInv"));
            }
        }
        return params;
    }

    public List<Parameter> getRelativeRateParameters() {
        List<Parameter> allMus = new ArrayList<>();
        int[] weights = getPartitionCodonWeights();
        if (getCodonPartitionCount() > 1) {
            if (codonHeteroPattern.equals("123")) {
                Parameter parameter = getParameter("CP1." + (options.useNuRelativeRates() ? "nu" : "mu"));
                parameter.setDimensionWeight(weights[0]);
                allMus.add(parameter);

                parameter = getParameter("CP2." + (options.useNuRelativeRates() ? "nu" : "mu"));
                parameter.setDimensionWeight(weights[1]);
                allMus.add(parameter);

                parameter = getParameter("CP3." + (options.useNuRelativeRates() ? "nu" : "mu"));
                parameter.setDimensionWeight(weights[2]);
                allMus.add(parameter);
            } else if (codonHeteroPattern.equals("112")) {
                Parameter parameter = getParameter("CP1+2." + (options.useNuRelativeRates() ? "nu" : "mu"));
                parameter.setDimensionWeight(weights[0]);
                allMus.add(parameter);

                parameter = getParameter("CP3." + (options.useNuRelativeRates() ? "nu" : "mu"));
                parameter.setDimensionWeight(weights[1]);
                allMus.add(parameter);
            } else {
                throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
            }
        } else {
            Parameter mu = getParameter(options.useNuRelativeRates() ? "nu" : "mu");
            mu.setDimensionWeight(weights[0]);
            allMus.add(mu);
        }
        return allMus;
    }

    private void addFrequencyParams(List<Parameter> params, boolean includeRelativeRates) {
        if (frequencyPolicy == FrequencyPolicyType.ESTIMATED) {
            if (includeRelativeRates && unlinkedSubstitutionModel && unlinkedFrequencyModel) {
                if (codonHeteroPattern.equals("123")) {
                    params.add(getParameter("CP1.frequencies"));
                    params.add(getParameter("CP2.frequencies"));
                    params.add(getParameter("CP3.frequencies"));
                } else if (codonHeteroPattern.equals("112")) {
                    params.add(getParameter("CP1+2.frequencies"));
                    params.add(getParameter("CP3.frequencies"));
                } else {
                    throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                }
            } else {
                params.add(getParameter("frequencies"));
            }

        }
    }

    @Override
    public List<Operator> selectOperators(List<Operator> operators) {

        List<Operator> ops = new ArrayList<Operator>();

        switch (getDataType().getType()) {
            case DataType.NUCLEOTIDES:

                if (hasCodonPartitions() && unlinkedSubstitutionModel) {
                    if (codonHeteroPattern.equals("123")) {
                        switch (nucSubstitutionModel) {
                            case JC:
                                break;
                            case HKY:
                                ops.add(getOperator("CP1.kappa"));
                                ops.add(getOperator("CP2.kappa"));
                                ops.add(getOperator("CP3.kappa"));
                                break;

                            case TN93:
                                ops.add(getOperator("CP1.kappa1"));
                                ops.add(getOperator("CP2.kappa1"));
                                ops.add(getOperator("CP3.kappa1"));
                                ops.add(getOperator("CP1.kappa2"));
                                ops.add(getOperator("CP2.kappa2"));
                                ops.add(getOperator("CP3.kappa2"));
                                break;

                            case GTR:

                                for (int i = 1; i <= 3; i++) {
                                    if (options.useNewGTR()) {
                                        ops.add(getOperator("CP" + i + ".deltaGTR"));
                                    } else {
                                        for (String rateName : GTR_RATE_NAMES) {
                                            ops.add(getOperator("CP" + i + "." + rateName));
                                        }
                                    }
                                }
                                break;

                            default:
                                throw new IllegalArgumentException("Unknown nucleotides substitution model");
                        }

                    } else if (codonHeteroPattern.equals("112")) {
                        switch (nucSubstitutionModel) {
                            case JC:
                                break;
                            case HKY:
                                ops.add(getOperator("CP1+2.kappa"));
                                ops.add(getOperator("CP3.kappa"));
                                break;

                            case TN93:
                                ops.add(getOperator("CP1+2.kappa1"));
                                ops.add(getOperator("CP3.kappa1"));
                                ops.add(getOperator("CP1+2.kappa2"));
                                ops.add(getOperator("CP3.kappa2"));
                                break;

                            case GTR:
                                if (options.useNewGTR()) {
                                    ops.add(getOperator("CP1+2.deltaGTR"));
                                    ops.add(getOperator("CP3.deltaGTR"));
                                } else {
                                    for (String rateName : GTR_RATE_NAMES) {
                                        ops.add(getOperator("CP1+2." + rateName));
                                    }
                                    for (String rateName : GTR_RATE_NAMES) {
                                        ops.add(getOperator("CP3." + rateName));
                                    }
                                }
                                break;

                            default:
                                throw new IllegalArgumentException("Unknown nucleotides substitution model");
                        }

                    } else {
                        throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                    }

                } else { // no codon partitioning, or unlinkedSubstitutionModel
                    switch (nucSubstitutionModel) {
                        case JC:
                            break;
                        case HKY:
                            ops.add(getOperator("kappa"));
                            break;

                        case TN93:
                            ops.add(getOperator("kappa1"));
                            ops.add(getOperator("kappa2"));
                            break;

                        case GTR:
                            if (options.useNewGTR()) {
                                ops.add(getOperator("deltaGTR"));
                            } else {
                                for (String rateName : GTR_RATE_NAMES) {
                                    ops.add(getOperator(rateName));
                                }
                            }
                            break;

                        default:
                            throw new IllegalArgumentException("Unknown nucleotides substitution model");
                    }
                }

                // only AMINO_ACIDS not addFrequency
                addFrequencyOps(ops);
                break;

            case DataType.AMINO_ACIDS:
                break;

            case DataType.TWO_STATES:
            case DataType.COVARION:
                switch (binarySubstitutionModel) {
                    case BIN_SIMPLE:
                    case BIN_DOLLO:
                        break;

                    case BIN_COVARION:
                        ops.add(getOperator("bcov.alpha"));
                        ops.add(getOperator("bcov.s"));
                        ops.add(getOperator("hfrequencies"));
                        break;

                    default:
                        throw new IllegalArgumentException("Unknown binary substitution model");
                }

                // only AMINO_ACIDS not addFrequency
                addFrequencyOps(ops);
                break;

            case DataType.GENERAL:
                break;

            case DataType.CONTINUOUS:
                break;

            case DataType.DUMMY:
                //Do nothing
                break;
            case DataType.TREE:
                //Do nothing
                // System.out.println("Tree data type ParitionSubstutionModel selectOperators");
                break;

            default:
                throw new IllegalArgumentException("Unknown data type");
        }

        // if gamma do shape move
        if (gammaHetero) {
            if (hasCodonPartitions() && unlinkedHeterogeneityModel) {
                if (codonHeteroPattern.equals("123")) {
                    ops.add(getOperator("CP1.alpha"));
                    ops.add(getOperator("CP2.alpha"));
                    ops.add(getOperator("CP3.alpha"));
                } else if (codonHeteroPattern.equals("112")) {
                    ops.add(getOperator("CP1+2.alpha"));
                    ops.add(getOperator("CP3.alpha"));
                } else {
                    throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                }
            } else {
                ops.add(getOperator("alpha"));
            }
        }
        // if pinv do pinv move
        if (invarHetero) {
            String name = (options.usePInvRandomWalk() ? "rwPInv" : "uniformPInv");
            if (hasCodonPartitions() && unlinkedHeterogeneityModel) {
                if (codonHeteroPattern.equals("123")) {
                    ops.add(getOperator("CP1." + name));
                    ops.add(getOperator("CP2." + name));
                    ops.add(getOperator("CP3." + name));
                } else if (codonHeteroPattern.equals("112")) {
                    ops.add(getOperator("CP1+2." + name));
                    ops.add(getOperator("CP3." + name));
                } else {
                    throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                }
            } else {
                ops.add(getOperator(name));
            }
        }

        if (options.operatorSetType != OperatorSetType.CUSTOM) {
            // unless a custom mix has been chosen these operators should be off if AMTK is on
            for (Operator op : ops) {
                if (op.getParameter1().isAdaptiveMultivariateCompatible) {
                    op.setUsed(options.operatorSetType != OperatorSetType.ADAPTIVE_MULTIVARIATE);
                }
            }
        }

        operators.addAll(ops);

        return operators;
    }

    private void addFrequencyOps(List<Operator> ops) {
        if (frequencyPolicy == FrequencyPolicyType.ESTIMATED) {
            if (hasCodonPartitions() && unlinkedSubstitutionModel && unlinkedFrequencyModel) {
                if (codonHeteroPattern.equals("123")) {
                    ops.add(getOperator("CP1.frequencies"));
                    ops.add(getOperator("CP2.frequencies"));
                    ops.add(getOperator("CP3.frequencies"));
                } else if (codonHeteroPattern.equals("112")) {
                    ops.add(getOperator("CP1+2.frequencies"));
                    ops.add(getOperator("CP3.frequencies"));
                } else {
                    throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                }
            } else {
                ops.add(getOperator("frequencies"));
            }
        }
    }

    /**
     * @return true either if the options have more than one partition or any partition is
     * broken into codon positions.
     */
    public boolean hasCodonPartitions() {
        return getCodonPartitionCount() > 1;
    }

    public int getCodonPartitionCount() {
        if (codonHeteroPattern == null || codonHeteroPattern.equals("111")) {
            return 1;
        }
        if (codonHeteroPattern.equals("123")) {
            return 3;
        }
        if (codonHeteroPattern.equals("112")) {
            return 2;
        }
        throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
    }

    public void addWeightsForPartition(AbstractPartitionData partition, int[] weights, int offset) {
        int n = partition.getSiteCount();
        int codonCount = n / 3;
        int remainder = n % 3;
        if (codonHeteroPattern == null || codonHeteroPattern.equals("111")) {
            weights[offset] += n;
            return;
        }
        if (codonHeteroPattern.equals("123")) {
            weights[offset] += codonCount + (remainder > 0 ? 1 : 0);
            weights[offset + 1] += codonCount + (remainder > 1 ? 1 : 0);
            weights[offset + 2] += codonCount;
            return;
        }
        if (codonHeteroPattern.equals("112")) {
            weights[offset] += codonCount * 2 + remainder; // positions 1 + 2
            weights[offset + 1] += codonCount; // position 3
            return;
        }
        throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
    }

    /**
     * This returns an integer vector of the number of sites in each partition (including any codon partitions). These
     * are strictly in the same order as the 'mu' relative rates are listed.
     *
     * @return weights for each partition model
     */
    public int[] getPartitionCodonWeights() {
        int[] weights = new int[getCodonPartitionCount()];

        int k = 0;
        for (AbstractPartitionData partition : options.getDataPartitions(this)) {
            if (partition.getPartitionSubstitutionModel() == this) {
                addWeightsForPartition(partition, weights, k);
            }
        }
        k += getCodonPartitionCount();

        assert (k == weights.length);

        return weights;
    }

    ///////////////////////////////////////////////////////

    public NucModelType getNucSubstitutionModel() {
        return nucSubstitutionModel;
    }

    public void setNucSubstitutionModel(NucModelType nucSubstitutionModel) {
        this.nucSubstitutionModel = nucSubstitutionModel;
    }

    public AminoAcidModelType getAaSubstitutionModel() {
        return aaSubstitutionModel;
    }

    public void setAaSubstitutionModel(AminoAcidModelType aaSubstitutionModel) {
        this.aaSubstitutionModel = aaSubstitutionModel;
    }

    public BinaryModelType getBinarySubstitutionModel() {
        return binarySubstitutionModel;
    }

    public void setBinarySubstitutionModel(BinaryModelType binarySubstitutionModel) {
        this.binarySubstitutionModel = binarySubstitutionModel;
    }

    public DiscreteSubstModelType getDiscreteSubstType() {
        return discreteSubstType;
    }

    public void setDiscreteSubstType(DiscreteSubstModelType discreteSubstType) {
        this.discreteSubstType = discreteSubstType;
    }

    public ContinuousSubstModelType getContinuousSubstModelType() {
        return continuousSubstModelType;
    }

    public ContinuousModelExtensionType getContinuousExtensionType() {
        return continuousExtensionType;
    }

    public void setContinuousSubstModelType(final ContinuousSubstModelType continuousSubstModelType) {
        this.continuousSubstModelType = continuousSubstModelType;
    }

    public void setContinuousExtensionType(final ContinuousModelExtensionType extensionType) {
        this.continuousExtensionType = extensionType;
    }

    public void setContinuousLatentDimension(final int k) {
        this.continuousTraitCount = k;
    }

    public void setIsLatitudeLongitude(boolean latitudeLongitude) {
        isLatitudeLongitude = latitudeLongitude;
    }

    public boolean isLatitudeLongitude() {
        return isLatitudeLongitude;
    }

    public void setIsIndependent(boolean isIndependent) {
        this.isIndependent = isIndependent;
    }

    public boolean isIndependent() {
        return isIndependent;
    }

    public void setJitterWindow(double jitterWindow) {
        this.jitterWindow = jitterWindow;
    }

    public double getJitterWindow() {
        return jitterWindow;
    }

    public int getContinuousTraitCount() {
        return continuousTraitCount;
    }

    public int getContinuousTraitDimension() {
        if (isIndependent) {
            return 1;
        }

        return getContinuousTraitCount();
    }

    public int getExtendedTraitCount() {
        return extendedTraitCount;
    }

    public boolean isActivateBSSVS() {
        return discreteSubstType != DiscreteSubstModelType.GLM_SUBST && activateBSSVS;
    }

    public void setActivateBSSVS(boolean activateBSSVS) {
        this.activateBSSVS = activateBSSVS;
    }

    public FrequencyPolicyType getFrequencyPolicy() {
        return frequencyPolicy;
    }

    public void setFrequencyPolicy(FrequencyPolicyType frequencyPolicy) {
        this.frequencyPolicy = frequencyPolicy;
    }

    public boolean isGammaHetero() {
        return gammaHetero;
    }

    public void setGammaHetero(boolean gammaHetero) {
        this.gammaHetero = gammaHetero;
    }

    public int getGammaCategories() {
        return gammaCategories;
    }

    public void setGammaCategories(int gammaCategories) {
        this.gammaCategories = gammaCategories;
    }

    public boolean isInvarHetero() {
        return invarHetero;
    }

    public void setInvarHetero(boolean invarHetero) {
        this.invarHetero = invarHetero;
    }

    public boolean isGammaHeteroEqualWeights() {
        return equalWeights;
    }

    public void setGammaHeteroEqualWeights(boolean equalWeights) {
        this.equalWeights = equalWeights;
    }

    public String getCodonHeteroPattern() {
        return codonHeteroPattern;
    }

    public void setCodonHeteroPattern(String codonHeteroPattern) {
        this.codonHeteroPattern = codonHeteroPattern;
    }

    /**
     * @return true if the rate matrix parameters are unlinked across codon positions
     */
    public boolean isUnlinkedSubstitutionModel() {
        return unlinkedSubstitutionModel;
    }

    public void setUnlinkedSubstitutionModel(boolean unlinkedSubstitutionModel) {
        this.unlinkedSubstitutionModel = unlinkedSubstitutionModel;
    }

    public boolean isUnlinkedHeterogeneityModel() {
        return unlinkedHeterogeneityModel;
    }

    public void setUnlinkedHeterogeneityModel(boolean unlinkedHeterogeneityModel) {
        this.unlinkedHeterogeneityModel = unlinkedHeterogeneityModel;
    }

    public boolean isUnlinkedFrequencyModel() {
        return unlinkedFrequencyModel;
    }

    public void setUnlinkedFrequencyModel(boolean unlinkedFrequencyModel) {
        this.unlinkedFrequencyModel = unlinkedFrequencyModel;
    }

    public boolean isDolloModel() {
        return dolloModel;
    }

    public void setDolloModel(boolean dolloModel) {
        this.dolloModel = dolloModel;
    }

    public boolean isUseAmbiguitiesTreeLikelihood() {
        return useAmbiguitiesTreeLikelihood;
    }

    public void setUseAmbiguitiesTreeLikelihood(boolean useAmbiguitiesTreeLikelihood) {
        this.useAmbiguitiesTreeLikelihood = useAmbiguitiesTreeLikelihood;
    }

    public String getPrefix() {
        String prefix = "";
        if (options.getPartitionSubstitutionModels(Nucleotides.INSTANCE).size() +
                options.getPartitionSubstitutionModels(AminoAcids.INSTANCE).size() > 1) {
            // There is more than one active partition model, or doing species analysis
            prefix += getName() + ".";
        }
        return prefix;
    }

    public String getPrefix(DataType dataType) {
        String prefix = "";
        if (options.getPartitionSubstitutionModels(dataType).size() > 1) {
            // There is more than one active partition model, or doing species analysis
            prefix += getName() + ".";
        }
        return prefix;
    }

    public String getPrefix(int codonPartitionNumber) {
        String prefix = "";
        prefix += getPrefix();
        prefix += getPrefixCodon(codonPartitionNumber);
        return prefix;
    }

    public String getPrefixCodon(int codonPartitionNumber) {
        String prefix = "";
        if (getCodonPartitionCount() > 1 && codonPartitionNumber > 0) {
            if (getCodonHeteroPattern().equals("123")) {
                prefix += "CP" + codonPartitionNumber + ".";
            } else if (getCodonHeteroPattern().equals("112")) {
                if (codonPartitionNumber == 1) {
                    prefix += "CP1+2.";
                } else {
                    prefix += "CP3.";
                }
            } else {
                throw new IllegalArgumentException("unsupported codon hetero pattern");
            }

        }
        return prefix;
    }

    /**
     * returns the union of the set of states for all traits using this discrete CTMC model
     *
     * @return
     */
    public Set<String> getDiscreteStateSet() {
        Set<String> states = new HashSet<String>();
        for (AbstractPartitionData partition : options.getDataPartitions(this)) {
            if (partition.getTraits() != null) {
                states.addAll(partition.getTraits().get(0).getStatesOfTrait(options.taxonList));
            }
        }
        return states;
    }

    public void copyFrom(PartitionSubstitutionModel source) {
        nucSubstitutionModel = source.nucSubstitutionModel;
        aaSubstitutionModel = source.aaSubstitutionModel;
        binarySubstitutionModel = source.binarySubstitutionModel;
        discreteSubstType = source.discreteSubstType;
        continuousSubstModelType = source.continuousSubstModelType;

        activateBSSVS = source.activateBSSVS;
        useAmbiguitiesTreeLikelihood = source.useAmbiguitiesTreeLikelihood;

        frequencyPolicy = source.frequencyPolicy;
        gammaHetero = source.gammaHetero;
        gammaCategories = source.gammaCategories;
        equalWeights = source.equalWeights;
        invarHetero = source.invarHetero;
        codonHeteroPattern = source.codonHeteroPattern;
        unlinkedSubstitutionModel = source.unlinkedSubstitutionModel;
        unlinkedHeterogeneityModel = source.unlinkedHeterogeneityModel;
        unlinkedFrequencyModel = source.unlinkedFrequencyModel;

        dolloModel = source.dolloModel;
    }

    @Override
    public String toString() {
        return getName();
    }

}