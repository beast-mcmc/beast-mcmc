/*
 * BeautiOptions.java
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

import java.util.*;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.components.ancestralstates.AncestralStatesComponentOptions;
import dr.app.beauti.components.continuous.ContinuousComponentOptions;
import dr.app.beauti.components.discrete.DiscreteTraitsComponentOptions;
import dr.app.beauti.mcmcpanel.MCMCPanel;
import dr.app.beauti.types.OperatorSetType;
import dr.app.beauti.types.TreePriorType;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.DataType;
import dr.evolution.util.Date;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evoxml.util.DateUnitsType;
import dr.inference.operators.OperatorSchedule;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class BeautiOptions extends ModelOptions {

    // Switches to a dirichlet prior & delta exchange on nus
    private static final boolean NEW_RELATIVE_RATE_PARAMETERIZATION = true;

    // Switches to a dirichlet prior & delta exchange on relative rates
    private static final boolean NEW_GTR_PARAMETERIZATION = true;

    // Makes sets the initial state of PartitionClockModel.continuousQuantile
    public static final boolean DEFAULT_QUANTILE_RELAXED_CLOCK = false;

    // Uses a logit transformed random walk on pInv parameters
    private static final boolean LOGIT_PINV_KERNEL = true;

    // Switches to a dirichlet prior & delta exchange on state frequencies
    private static final boolean FREQUENCIES_DIRICHLET_PRIOR = true;

    // Switches to a gamma prior on population size parameters
    private static final boolean POPSIZE_GAMMA_PRIOR = true;

    private final String DEFAULT_NAME = "default";

    private static final long serialVersionUID = -3676802825545741012L;

    public BeautiOptions() {
        this(new ComponentFactory[]{});
    }

    public BeautiOptions(ComponentFactory[] components) {

        // Install all the component's options from the given list of factories:
        registerComponents(components);
    }

    /**
     * This will register the list of components if not already there...
     *
     * @param components
     */
    public void registerComponents(ComponentFactory[] components) {

        // Install all the component's options from the given list of factories:
        for (ComponentFactory component : components) {
            if (!hasComponent(component)) {
                addComponent(component.createOptions(this));
            }
        }
    }


    /**
     * resets the options to the initial conditions
     */
    public void reset() {
        // Data options
        taxonList = null;
        taxonSets.clear();
        taxonSetsMono.clear();
        taxonSetsIncludeStem.clear();
        taxonSetsTreeModel.clear();

        // meanDistance = 1.0;
        datesUnits = DateUnitsType.YEARS;
        datesDirection = DateUnitsType.FORWARDS;
        maximumTipHeight = 0.0;
        translation = 0;

        dataPartitions.clear();
        traits.clear();
        // partitionModels.clear();
        // partitionTreeModels.clear();
        // partitionTreePriors.clear();
        // partitionClockTreeLinks.clear();
        // activedSameTreePrior = null;
        // shareSameTreePrior = true;
        userTrees.clear();

        // rateOptionClockModel = FixRateType.FIX_FIRST_PARTITION;
        // meanSubstitutionRate = 1.0;
        unlinkPartitionRates = true;

        units = Units.Type.SUBSTITUTIONS;

        // Operator schedule options
        optimizationTransform = OperatorSchedule.DEFAULT_TRANSFORM;

        // MCMC options
        chainLength = 10000000;
        logEvery = 1000;
        echoEvery = 1000;
        checkpointEvery = 1000000;
        burnIn = 100000;
        fileName = null;
        autoOptimize = true;
        performTraceAnalysis = false;
        generateCSV = true;  // until/if a button
        samplePriorOnly = false;

        fileNameStem = MCMCPanel.DEFAULT_FILE_NAME_STEM;
        logFileName = null;
        checkpointFileName = null;
        allowOverwriteLog = false;
//        mapTreeLog = false;
//        mapTreeFileName = null;
        treeFileName.clear();
        substTreeLog = false;
        substTreeFileName.clear();
        operatorAnalysis = false;
        operatorAnalysisFileName = null;

        globalModelOptions = new GlobalModelOptions(this);
        clockModelOptions = new ClockModelOptions(this);
        treeModelOptions = new TreeModelOptions(this);
//        priorOptions = new PriorOptions(this);

//        traitsOptions = new TraitsOptions(this);

        parameters.clear();
        operators.clear();
        statistics.clear();

        clearDataPartitionCaches();
    }

    public void selectTaxonSetsStatistics(List<Parameter> params) {

        if (taxonSets != null) {
            for (Taxa taxa : taxonSets) {
                params.add(updateTMRCAStatistic(taxa));
            }
        } else {
            System.err.println("TaxonSets are null");
        }
    }

    public Parameter updateTMRCAStatistic(Taxa taxa) {
        Parameter statistic = statistics.get(taxa);
        PartitionTreeModel treeModel = taxonSetsTreeModel.get(taxa);
        if (statistic == null) {
            // default scaleType = PriorScaleType.NONE; priorType = PriorType.NONE_TREE_PRIOR
            statistic = new Parameter.Builder(taxa.getId(), "tmrca statistic for taxon set")
                    .isStatistic(true).isNodeHeight(true)
                    .partitionOptions(treeModel).initial(Double.NaN).isNonNegative(true).build();
            statistic.setPrefix(treeModel.getPrefix());
            statistic.taxaId = taxa.getId();

            statistics.put(taxa, statistic);

        } else {
            statistic.setOptions(treeModel); // keep consistent to taxonSetsTreeModel
            statistic.setPrefix(treeModel.getPrefix()); // keep prefix consistent after link/unlink tree
            PartitionTreePrior treePrior = treeModel.getPartitionTreePrior();
            statistic.isCalibratedYule = treePrior.getNodeHeightPrior() == TreePriorType.YULE_CALIBRATION
                    && taxonSetsMono.get(taxa);
            statistic.taxaId = taxa.getId();
        }
        return statistic;
    }

    public boolean renameTMRCAStatistic(Taxa taxonSet) {
        Parameter statistic = statistics.get(taxonSet);
        if (statistic != null) {
            PartitionTreeModel treeModel = taxonSetsTreeModel.get(taxonSet);
            statistic.taxaId = treeModel.getPrefix() + taxonSet.getId();
            return true;
        } else {
            return false;
        }
    }

    @Override
    public void initModelParametersAndOpererators() {

    }

    public List<Parameter> selectParameters() {
        return selectParameters(new ArrayList<Parameter>());
    }


    /**
     * return an list of parameters that are required
     *
     * @return the parameter list
     */
    @Override
    public List<Parameter> selectParameters(List<Parameter> parameters) {

        globalModelOptions.selectParameters(parameters);

        selectTaxonSetsStatistics(parameters); // have to be before clockModelOptions.selectParameters(parameters);

        for (PartitionSubstitutionModel model : getPartitionSubstitutionModels()) {
//          parameters.addAll(model.getParameters(multiplePartitions));
            model.selectParameters(parameters);
        }

        for (PartitionClockModel model : getPartitionClockModels()) {
            Set<PartitionSubstitutionModel> substitutionModels = new LinkedHashSet<>();
            for (AbstractPartitionData partition : getDataPartitions()) {
                if (partition.getPartitionClockModel() == model) {
                    PartitionSubstitutionModel substitutionModel = partition.getPartitionSubstitutionModel();
                    if(substitutionModel!=null){
                        substitutionModels.add(substitutionModel);
                    }
                    // Tree data doesn't have subsitution model at this point
                }
            }

            // collect all the relative rate paremeters (partition rates and codon position rates)
            ArrayList<Parameter> relativeRateParameters = new ArrayList<Parameter>();
            for (PartitionSubstitutionModel substitutionModel : substitutionModels) {
                relativeRateParameters.addAll(substitutionModel.getRelativeRateParameters());
            }
            Parameter allMuNus = model.getParameter(useNuRelativeRates() ? "allNus" : "allMus");
            allMuNus.clearSubParameters();
            if (relativeRateParameters.size() > 1) {

                int totalWeight = 0;
                for (Parameter mu : relativeRateParameters) {
                    allMuNus.addSubParameter(mu);
                    totalWeight += mu.getDimensionWeight();
                }
                parameters.add(allMuNus);

                // add the total weight of all mus/nus to the allMus parameter
                allMuNus.setDimensionWeight(totalWeight);
                allMuNus.maintainedSum = (useNuRelativeRates() ? 1.0 : relativeRateParameters.size());
                allMuNus.isMaintainedSum = true;
            }

            model.selectParameters(parameters);
        }
        clockModelOptions.selectParameters(parameters);

        for (PartitionTreeModel tree : getPartitionTreeModels()) {
            tree.selectParameters(parameters);
        }
        treeModelOptions.selectParameters(parameters);

        for (PartitionTreePrior prior : getPartitionTreePriors()) {
            //refresh prior choices
            prior.alternatePopulationSizePriors();

            prior.selectParameters(parameters);
        }

//        for (TraitData trait : getTraitsList()) { // all traits including locations
//            if (!trait.getName().equalsIgnoreCase(TraitData.Traits.TRAIT_SPECIES.toString()))
//        	   trait.gets.selectParameters(parameters);
//        }

        selectComponentParameters(this, parameters);

        selectComponentStatistics(this, parameters);

//        priorOptions.selectParameters(parameters);

        return parameters;
    }

    /**
     * return an list of operators that are required
     *
     * @return the operator list
     */
    public List<Operator> selectOperators() {
        return selectOperators(new ArrayList<Operator>());
    }

    @Override
    public List<Operator> selectOperators(List<Operator> ops) {
        globalModelOptions.selectOperators(ops);

        for (PartitionSubstitutionModel model : getPartitionSubstitutionModels()) {
            model.selectOperators(ops);
        }

        for (PartitionClockModel model : getPartitionClockModels()) {
            model.selectOperators(ops);
        }
        clockModelOptions.selectOperators(ops);

        for (PartitionTreeModel tree : getPartitionTreeModels()) {
            tree.selectOperators(ops);
        }
        treeModelOptions.selectOperators(ops);

        for (PartitionTreePrior prior : getPartitionTreePriors()) {
            prior.selectOperators(ops);
        }

//        for (TraitData trait : getTraitsList()) { // all traits including locations
//        	if (!trait.getName().equalsIgnoreCase(TraitData.Traits.TRAIT_SPECIES.toString()))
//                trait.getTraitData().selectOperators(ops);
//        }

        selectComponentOperators(this, ops);

        // no remove operators for parameters that are part of a joint prior...
        List<Operator> toRemove = new ArrayList<Operator>();
        for (Operator operator : ops) {
            if ((operator.getParameter1() != null && operator.getParameter1().isLinked) ||
                    (operator.getParameter2() != null && operator.getParameter2().isLinked)) {
                toRemove.add(operator);
            }
        }

        ops.removeAll(toRemove);

        return ops;
    }

    @Override
    public String getPrefix() {
        return "";
    }

    public Operator getOperator(Parameter parameter) {
        for (Operator operator : selectOperators()) {
            if (operator.getParameter1() == parameter || operator.getParameter2() == parameter) {
                return operator;
            }
        }
        return null;
    }


    public boolean hasData() {
        return dataPartitions.size() > 0;
    }

    public boolean contains(DataType dataType) {
        for (AbstractPartitionData pd : dataPartitions) {
            if (pd.getDataType().getType() == dataType.getType()) {
                return true;
            }
        }
        return false;
    }

    public boolean hasPartitionData(String name) {
        for (AbstractPartitionData pd : dataPartitions) {
            if (name.equalsIgnoreCase(pd.getName())) {
                return true;
            }
        }
        return false;
    }

    public PartitionData getPartitionData(Alignment alignment) {
        for (PartitionData partition : getPartitionData()) {
            if (partition.getAlignment() == alignment)
                return partition;
        }
        return null;
    }

    /**
     * exclude traits
     */
    public List<PartitionData> getPartitionData() {
        List<PartitionData> pdList = new ArrayList<>();
        for (AbstractPartitionData partition : dataPartitions) {
            if (partition instanceof PartitionData && partition.getTraits() == null) {
                pdList.add((PartitionData) partition);
            }
        }
        return pdList;
    }

    public List<AbstractPartitionData> getDataPartitions() {
        return dataPartitions;
    }

    public List<AbstractPartitionData> getTraitPartitions(TraitData trait) {
        List<AbstractPartitionData> pdList = new ArrayList<AbstractPartitionData>();
        for (AbstractPartitionData pd : dataPartitions) {
            if (pd.getTraits() != null && pd.getTraits().contains(trait)) {
                pdList.add(pd);
            }
        }
        return pdList;
    }

    public List<AbstractPartitionData> getDataPartitions(DataType dataType) {
        List<AbstractPartitionData> pdList = new ArrayList<AbstractPartitionData>();
        for (AbstractPartitionData pd : dataPartitions) {
            if (pd.getDataType().getType() == dataType.getType()) {
                pdList.add(pd);
            }
        }
        return pdList;
    }

    public List<AbstractPartitionData> getDataPartitionsExceptTrees() {
        List<AbstractPartitionData> pdList = new ArrayList<AbstractPartitionData>();
        for (AbstractPartitionData pd : dataPartitions) {
            if (pd.getDataType().getType() != DataType.TREE) {
                pdList.add(pd);
            }
        }
        return pdList;
    }

    public List<AbstractPartitionData> getDataPartitions(PartitionOptions model) {
        if (model instanceof PartitionSubstitutionModel) {
            return getDataPartitions((PartitionSubstitutionModel) model);
        } else if (model instanceof PartitionClockModel) {
            return getDataPartitions((PartitionClockModel) model);
        } else if (model instanceof PartitionTreeModel) {
            return getDataPartitions((PartitionTreeModel) model);
        } else if (model instanceof PartitionTreePrior) {
            return getDataPartitions((PartitionTreePrior) model);
        } else {
            return null;
        }
    }

    private final Map<PartitionSubstitutionModel, List<AbstractPartitionData>> psmCache = new HashMap<PartitionSubstitutionModel, List<AbstractPartitionData>>();

    public List<AbstractPartitionData> getDataPartitions(PartitionSubstitutionModel model) {
        List<AbstractPartitionData> pdList = psmCache.get(model);

        if (pdList == null) {
            pdList = new ArrayList<AbstractPartitionData>();

            for (AbstractPartitionData pd : dataPartitions) {
                if (pd.getPartitionSubstitutionModel() == model) {
                    pdList.add(pd);
                }
            }

            psmCache.put(model, pdList);
        }
        return pdList;
    }

    private final Map<PartitionClockModel, List<AbstractPartitionData>> pcmCache = new HashMap<PartitionClockModel, List<AbstractPartitionData>>();
    private final Map<PartitionTreeModel, List<AbstractPartitionData>> ptmCache = new HashMap<PartitionTreeModel, List<AbstractPartitionData>>();

    public List<AbstractPartitionData> getDataPartitions(PartitionTreeModel model) {
        List<AbstractPartitionData> pdList = ptmCache.get(model);

        if (pdList == null) {
            pdList = new ArrayList<AbstractPartitionData>();

            for (AbstractPartitionData pd : dataPartitions) {
                if (pd.getPartitionTreeModel() == model) {
                    pdList.add(pd);
                }
            }

            ptmCache.put(model, pdList);
        }
        return pdList;
    }

    private final Map<PartitionTreePrior, List<AbstractPartitionData>> ptpCache = new HashMap<PartitionTreePrior, List<AbstractPartitionData>>();

    public List<AbstractPartitionData> getDataPartitions(PartitionTreePrior prior) {
        List<AbstractPartitionData> pdList = ptpCache.get(prior);

        if (pdList == null) {
            pdList = new ArrayList<AbstractPartitionData>();

            for (AbstractPartitionData pd : dataPartitions) {
                if (pd.getPartitionTreeModel().getPartitionTreePrior() == prior) {
                    pdList.add(pd);
                }
            }

            ptpCache.put(prior, pdList);
        }
        return pdList;
    }

    public List<AbstractPartitionData> getDataPartitions(PartitionClockModel clockModel) {
        List<AbstractPartitionData> pdList = pcmCache.get(clockModel);

        if (pdList == null) {
            pdList = new ArrayList<AbstractPartitionData>();

            for (AbstractPartitionData pd : dataPartitions) {
                if (pd.getPartitionClockModel() != null && pd.getPartitionClockModel() == clockModel) {
                    pdList.add(pd);
                }
            }

            pcmCache.put(clockModel, pdList);
        }
        return pdList;
    }

    public void clearDataPartitionCaches() {
        psmCache.clear();
        pcmCache.clear();
        ptmCache.clear();
        ptpCache.clear();
        psmlCache.clear();
        ptmlCache.clear();
        pcmlCache.clear();
    }

    public List<PartitionSubstitutionModel> getPartitionSubstitutionModels(DataType dataType) {
        List<PartitionSubstitutionModel> models = new ArrayList<PartitionSubstitutionModel>();
        for (PartitionSubstitutionModel model : getPartitionSubstitutionModels(dataPartitions)) {
            if (model.getDataType().getType() == dataType.getType()) {
                models.add(model);
            }
        }
        return models;
    }

    private final Map<List<? extends AbstractPartitionData>, List<PartitionSubstitutionModel>> psmlCache = new HashMap<List<? extends AbstractPartitionData>, List<PartitionSubstitutionModel>>();

    public List<PartitionSubstitutionModel> getPartitionSubstitutionModels(List<? extends AbstractPartitionData> givenDataPartitions) {
        List<PartitionSubstitutionModel> psmList = psmlCache.get(givenDataPartitions);

        if (psmList == null) {
            Set<PartitionSubstitutionModel> activeModels = new LinkedHashSet<PartitionSubstitutionModel>();

            for (AbstractPartitionData partition : givenDataPartitions) {
                PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();
                if (model != null) {
                    activeModels.add(model);
                }
            }

            psmList = new ArrayList<PartitionSubstitutionModel>(activeModels);

            psmlCache.put(givenDataPartitions, psmList);
        }

        return psmList;
    }

    public List<PartitionSubstitutionModel> getPartitionSubstitutionModels() {
        return getPartitionSubstitutionModels(dataPartitions);
    }

//    public List<PartitionSubstitutionModel> getPartitionTraitsSubstitutionModels() {
//        return getPartitionSubstitutionModels(getTraitsList());
//    }
//
//    public List<PartitionSubstitutionModel> getPartitionNonTraitsSubstitutionModels() {
//        return getPartitionSubstitutionModels(getNonTraitsDataList());
//    }


    // ++++++++++++++ Partition Clock Model ++++++++++++++

    private final Map<List<? extends AbstractPartitionData>, List<PartitionClockModel>> pcmlCache = new HashMap<List<? extends AbstractPartitionData>, List<PartitionClockModel>>();

    public List<PartitionClockModel> getPartitionClockModels(List<? extends AbstractPartitionData> givenDataPartitions) {
        List<PartitionClockModel> pcmList = pcmlCache.get(givenDataPartitions);

        if (pcmList == null) {
            Set<PartitionClockModel> activeModels = new LinkedHashSet<PartitionClockModel>();

            for (AbstractPartitionData partition : givenDataPartitions) {
                PartitionClockModel model = partition.getPartitionClockModel();
                if (model != null && !(partition.getDataType().getType() == DataType.CONTINUOUS)) {
                    activeModels.add(model);
                }
            }

            pcmList = new ArrayList<PartitionClockModel>(activeModels);

            pcmlCache.put(givenDataPartitions, pcmList);
        }

        return pcmList;
    }

    public List<PartitionClockModel> getPartitionClockModels(DataType dataType) {
        List<PartitionClockModel> models = new ArrayList<PartitionClockModel>();
        for (PartitionClockModel model : getPartitionClockModels()) {
            if (model.getDataType().getType() == dataType.getType()) {
                models.add(model);
            }
        }
        return models;
    }

    public List<PartitionClockModel> getPartitionClockModels() {
        return getPartitionClockModels(dataPartitions);
    }

    // ++++++++++++++ Partition Tree Model ++++++++++++++
//    public void addPartitionTreeModel(PartitionTreeModel tree) {
//
//        if (!partitionTreeModels.contains(tree)) {
//            partitionTreeModels.add(tree);
//        }
//    }

//    public List<PartitionTreeModel> getPartitionTreeModels() {
//        return partitionTreeModels;
//    }

    private final Map<List<? extends AbstractPartitionData>, List<PartitionTreeModel>> ptmlCache = new HashMap<List<? extends AbstractPartitionData>, List<PartitionTreeModel>>();

    public List<PartitionTreeModel> getPartitionTreeModels(List<? extends AbstractPartitionData> givenDataPartitions) {
        List<PartitionTreeModel> ptmList = ptmlCache.get(givenDataPartitions);

        if (ptmList == null) {
            Set<PartitionTreeModel> activeTrees = new LinkedHashSet<>();

            for (AbstractPartitionData partition : givenDataPartitions) {
                if (partition.getPartitionTreeModel() != null) {
                    activeTrees.add(partition.getPartitionTreeModel());
                }
            }

            ptmList = new ArrayList<>(activeTrees);

            ptmlCache.put(givenDataPartitions, ptmList);
        }
        return ptmList;
    }

    public List<PartitionTreeModel> getPartitionTreeModels() {
        return getPartitionTreeModels(dataPartitions);
    }

//    public List<PartitionTreeModel> getNonTraitPartitionTreeModels() {
//        return getPartitionTreeModels(getNonTraitsDataList());
//    }


    // ++++++++++++++ Partition Tree Prior ++++++++++++++

    public List<PartitionTreePrior> getPartitionTreePriors() {

        List<PartitionTreePrior> activeTrees = new ArrayList<PartitionTreePrior>();

        // # tree prior = 1 or # tree model
        for (PartitionTreeModel model : getPartitionTreeModels()) {
            PartitionTreePrior prior = model.getPartitionTreePrior();
            if (prior != null && (!activeTrees.contains(prior))) {
                activeTrees.add(prior);
            }
        }

        return activeTrees;
    }

    public List<PartitionTreeModel> getPartitionTreeModels(PartitionTreePrior treePrior) {

        List<PartitionTreeModel> activeTrees = new ArrayList<PartitionTreeModel>();

        for (PartitionTreeModel model : getPartitionTreeModels()) {
            PartitionTreePrior prior = model.getPartitionTreePrior();
            if (prior != null && model.getPartitionTreePrior() == treePrior && (!activeTrees.contains(model))) {
                activeTrees.add(model);
            }
        }

        return activeTrees;
    }

    /**
     * unlink all and copy the tree prior in selectedTreeModel to others
     * currently, tree prior name cannot be changed, but tree model name can, so that we have to use instance
     *
     * @param selectedTreeModel the selected tree model whose tree prior copied to others
     */
    public void unLinkTreePriors(PartitionTreeModel selectedTreeModel) {
        for (PartitionTreeModel model : getPartitionTreeModels()) {
            // because # tree prior = 1 or # tree model, prior here will be a same instance through all tree models
            PartitionTreePrior prior = model.getPartitionTreePrior();
            if (model == selectedTreeModel) {
                prior.setName(model.getName()); // keep name same as its tree model
            } else {
                PartitionTreePrior ptp = new PartitionTreePrior(this, model.getName(), prior);
                model.setPartitionTreePrior(ptp);
            }
        }
        clearDataPartitionCaches();
    }

    // link all to given treePrior
    public void linkTreePriors(PartitionTreePrior treePrior) {
        if (treePrior == null) treePrior = new PartitionTreePrior(this, getPartitionTreeModels().get(0));
        for (PartitionTreeModel model : getPartitionTreeModels()) {
            model.setPartitionTreePrior(treePrior);
        }
        clearDataPartitionCaches();
    }

    public boolean isShareSameTreePrior() {
        return getPartitionTreePriors().size() <= 1;
    }

    // ++++++++++++++ Partition Clock Model ++++++++++++++

    public double getAveWeightedMeanDistance(List<AbstractPartitionData> partitions) {
        double meanDistance = 0;
        double totalSite = 0;
        for (AbstractPartitionData partition : partitions) {
            meanDistance = meanDistance + partition.getMeanDistance() * partition.getSiteCount();
            totalSite = totalSite + partition.getSiteCount();
        }

        if (totalSite == 0) {
            return 0;
        } else {
            return meanDistance / totalSite;
        }
    }

    /**
     * check if all taxa are same across partitions.
     */
    public boolean hasIdenticalTaxa(List<AbstractPartitionData> partitionDataList) {
        TaxonList taxa = null;
        for (AbstractPartitionData partition : partitionDataList) {
            if (taxa == null) {
                taxa = partition.getTaxonList();
            } else {
                final TaxonList taxa1 = partition.getTaxonList();
                if (taxa1 != null) {
                    // Trait partitions are assumed to have the same taxon list as the
                    // tree they are using and so will return null.
                    if (taxa1.getTaxonCount() != taxa.getTaxonCount()) {
                        return false;
                    }
                    for (int k = 0; k < taxa1.getTaxonCount(); ++k) {
                        if (taxa.getTaxonIndex(taxa1.getTaxonId(k)) == -1) {
                            return false;
                        }
                    }
                }
            }
        }
        return true;
    }

    public boolean hasIdenticalTaxa() {
        return hasIdenticalTaxa(dataPartitions);
    }

    public void updateTaxonList() {
        taxonList.removeAllTaxa();
        for (AbstractPartitionData partition : dataPartitions) {
            if (!partition.isCreatedFromTrait()) {
                for (Taxon taxon : partition.getTaxonList()) {
                    if (!taxonList.contains(taxon)) {
                        taxonList.addTaxon(taxon);
                    }
                }
            }
        }
    }

    /**
     * given a list of BEAUti AbstractPartitionData, take the union of taxon list from them
     * but if partition instanceof PartitionPattern and taxon is masked in Patterns class, then not count.
     *
     * @param partitionDataList can be BEAUti PartitionData or PartitionPattern or both
     * @return num of taxon
     */
    public int getTaxonCount(List<AbstractPartitionData> partitionDataList) {
        if (partitionDataList == null || partitionDataList.isEmpty()) return 0;

        List<String> taxonNameList = new ArrayList<>();
        for (AbstractPartitionData partition : partitionDataList) {
            if (partition.getTaxonList() != null) { // not a trait partition
                for (Taxon t : partition.getTaxonList()) {
                    if (!taxonNameList.contains(t.getId())) {
                        taxonNameList.add(t.getId());
                    }
                }
            }
        }
        return taxonNameList.size();
    }

    // +++++++++++++ Traits +++++++++++++

//    public List<PartitionData> getNonTraitsDataList() {
//        List<PartitionData> nonTraitsData = new ArrayList<PartitionData>();
//        for (PartitionData partition : dataPartitions) {
//            if (partition.getTraitType() == null) {
//                nonTraitsData.add(partition);
//            }
//        }
//        return nonTraitsData;
//    }
//
//    public List<TraitData> getTraitsList() {
//        List<TraitData> traits = new ArrayList<TraitData>();
//        for (PartitionData partition : dataPartitions) {
//            if (partition.getTraitType() != null) {
//                traits.add((TraitData) partition);
//            }
//        }
//        return traits;
//    }
//
//    public List<TraitData> getDiscreteIntegerTraits() {
//        List<TraitData> traits = new ArrayList<TraitData>();
//        for (PartitionData partition : dataPartitions) {
//            if (partition.getTraitType() != null && partition.getTraitType() != TraitData.TraitType.CONTINUOUS) {
//                traits.add((TraitData) partition);
//            }
//        }
//        return traits;
//    }
//
//    public boolean hasDiscreteIntegerTraitsExcludeSpecies() { // exclude species at moment
//        return getDiscreteIntegerTraits().size() > 1
//                || (getDiscreteIntegerTraits().size() > 0 && (!traitExists(TraitData.TRAIT_SPECIES)));
//    }

    public int getIndexOfTrait(String traitName) {
        int i = 0;
        for (TraitData trait : traits) {
            if (trait.getName().equalsIgnoreCase(traitName)) {
                return i;
            }
            i++;
        }
        return -1;
    }

    public boolean traitExists(String traitName) {
        return getIndexOfTrait(traitName) != -1;
    }

    public int addTrait(TraitData newTrait) {
        int selRow = getIndexOfTrait(newTrait.getName());

        if (selRow == -1) {
            traits.add(newTrait);
            selRow = traits.size() - 1; // start 0
        }

        return selRow; // only for trait panel
    }

    public int createPartitionForTraits(String name, TraitData trait) {
        List<TraitData> traits = new ArrayList<TraitData>();
        traits.add(trait);
        return createPartitionForTraits(name, traits);
    }

    public int createPartitionForTraits(String name, List<TraitData> traits) {
        int selRow = -1;

        PartitionData partition = new PartitionData(this, name, traits.get(0).getFileName(), traits);
        dataPartitions.add(partition);
        selRow = dataPartitions.size() - 1;

        if (partition.getPartitionSubstitutionModel() == null) {
            PartitionSubstitutionModel substModel = new PartitionSubstitutionModel(this, partition.getName(),
                    partition);
            partition.setPartitionSubstitutionModel(substModel);
        }

        if (partition.getPartitionTreeModel() == null) {
            if (getPartitionTreeModels().isEmpty()) {
                PartitionTreeModel treeModel = new PartitionTreeModel(this, DEFAULT_NAME);
                partition.setPartitionTreeModel(treeModel);
                PartitionTreePrior ptp = new PartitionTreePrior(this, treeModel);
                treeModel.setPartitionTreePrior(ptp);
            } else {
                partition.setPartitionTreeModel(getPartitionTreeModels().get(0));// always use 1st tree
            }
        }

        if (partition.getPartitionClockModel() == null && partition.getDataType().getType() != DataType.CONTINUOUS) {
            // PartitionClockModel based on PartitionData
            PartitionClockModel pcm = new PartitionClockModel(this, partition.getName(), partition, partition.getPartitionTreeModel());
            partition.setPartitionClockModel(pcm);
        }

        updateTraitParameters(partition);

        return selRow; // only for trait panel
    }

    public boolean onlyContinuousPartitions() {
        if (dataPartitions.isEmpty()) {
            return false;
        }

        boolean onlyContinuous = true;

        for (AbstractPartitionData partition : dataPartitions) {
            if (partition.getDataType().getType() != DataType.CONTINUOUS) {
                onlyContinuous = false;
                break;
            }
        }

        return onlyContinuous;
    }

    private void updateTraitParameters(AbstractPartitionData partition) {
        if (partition.isCreatedFromTrait()) {
            ContinuousComponentOptions comp = (ContinuousComponentOptions) getComponentOptions(ContinuousComponentOptions.class);
            comp.createParameters(this);

            DiscreteTraitsComponentOptions comp2 = (DiscreteTraitsComponentOptions) getComponentOptions(DiscreteTraitsComponentOptions.class);
            comp2.createParameters(this);

            AncestralStatesComponentOptions comp3 = (AncestralStatesComponentOptions) getComponentOptions(AncestralStatesComponentOptions.class);
            comp3.setReconstructAtNodes(partition, true);
            comp3.setReconstructAtMRCA(partition, false);
        }
    }
    public void setClockAndTree(AbstractPartitionData partition) {

        PartitionTreeModel treeModel;

        // use same tree model and same tree prior in beginning
        if (getPartitionTreeModels().isEmpty()) {
            // PartitionTreeModel based on PartitionData
            treeModel = new PartitionTreeModel(this, DEFAULT_NAME);
            partition.setPartitionTreeModel(treeModel);

            // PartitionTreePrior always based on PartitionTreeModel
            PartitionTreePrior ptp = new PartitionTreePrior(this, treeModel);
            treeModel.setPartitionTreePrior(ptp);
        } else { //if (options.getPartitionTreeModels() != null) {
//                        && options.getPartitionTreeModels().size() == 1) {
            treeModel = getPartitionTreeModels().get(0); // same tree model,
            partition.setPartitionTreeModel(treeModel); // if same tree model, therefore same prior
        }

        // use same clock model in beginning, have to create after partition.setPartitionTreeModel(ptm);
        if (getPartitionClockModels(partition.getDataType()).isEmpty()) {
            // PartitionClockModel based on PartitionData
            PartitionClockModel pcm = new PartitionClockModel(this, DEFAULT_NAME, partition, treeModel);
            partition.setPartitionClockModel(pcm);
        } else { //if (options.getPartitionClockModels() != null) {
//                        && options.getPartitionClockModels().size() == 1) {
            PartitionClockModel pcm = getPartitionClockModels(partition.getDataType()).get(0);
            partition.setPartitionClockModel(pcm);
        }
    }

    public void renamePartition(AbstractPartitionData partition, String newName) {
        partition.setName(newName);
        updateTraitParameters(partition);
    }

    public void removeTrait(String traitName) {
        if (traitExists(traitName)) {
            clearTraitValues(traitName); // Clear trait values
            traits.remove(getTrait(traitName));
        }
    }

    public void clearTraitValues(String traitName) {
        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            taxonList.getTaxon(i).setAttribute(traitName, "");
        }
    }

    public TraitData getTrait(String traitName) {
        for (TraitData trait : traits) {
            if (trait.getName().equalsIgnoreCase(traitName))
                return trait;
        }
        return null;
    }

//    public boolean hasAlignmentPartition() {
//        for (AbstractPartitionData partition : dataPartitions) {
//            if (partition instanceof PartitionData) {
//                if (((PartitionData) partition).getAlignment() != null) {
//                    return true;
//                }
//            }
//        }
//        return false;
//    }


    public boolean hasDiscreteTrait() {
        for (TraitData traitData : traits) {
            if (traitData.getTraitType() == TraitData.TraitType.DISCRETE) {
                return true;
            }
        }
        return false;
    }

    public boolean hasDiscreteTraitPartition() {
        for (AbstractPartitionData partition : dataPartitions) {
            if (partition.getTraits() != null && partition.getTraits().get(0).getTraitType() == TraitData.TraitType.DISCRETE) {
                return true;
            }
        }
        return false;
    }

    public boolean hasContinuousTrait() {
        for (TraitData traitData : traits) {
            if (traitData.getTraitType() == TraitData.TraitType.CONTINUOUS) {
                return true;
            }
        }
        return false;
    }

    public boolean hasContinuousTraitPartition() {
        for (AbstractPartitionData partition : dataPartitions) {
            if (partition.getTraits() != null && partition.getTraits().get(0).getTraitType() == TraitData.TraitType.CONTINUOUS) {
                return true;
            }
        }
        return false;
    }

    public Set<String> getStatesForDiscreteModel(PartitionSubstitutionModel model) {
        Set<String> states = new TreeSet<String>();
        for (AbstractPartitionData partition : getDataPartitions(model)) {
            Set<String> newStates = partition.getTraits().get(0).getStatesOfTrait(taxonList);

            if (states.size() > 0) {
                Set<String> shared = new HashSet<String>(states);
                shared.retainAll(newStates);
                if (shared.size() == 0) {
                    throw new IllegalArgumentException("For discrete trait partitions to have a linked model they must share states");
                }
            }

            states.addAll(newStates);
        }

        if (states.size() < 1) throw new IllegalArgumentException("The number of states must be greater than 1");

        return states;
    }

    public static TraitData.TraitType guessTraitType(TaxonList taxa, String name) {
        TraitData.TraitType type = TraitData.TraitType.DISCRETE;
        return type;
    }

    // ++++++++++++++++++++ tree partition +++++++++++++++++

    public int createPartitionForTree(TreeHolder trees, String partitionName) {
        int selRow = -1;

        TreePartitionData partition = new TreePartitionData(this, partitionName, trees.getFileName(), trees);
        dataPartitions.add(partition);
        selRow = dataPartitions.size() - 1;



        if (partition.getPartitionTreeModel() == null) {
            PartitionTreeModel treeModel = new PartitionTreeModel(this, partition);
            PartitionTreePrior partitionTreePrior = new PartitionTreePrior(this, treeModel);
            treeModel.setPartitionTreePrior(partitionTreePrior);
            partition.setPartitionTreeModel(treeModel);// always use 1st tree
            treeModel.setEmpiricalTreesFilename(trees.getFileName());
        }
        if(partition.getPartitionClockModel() == null) {
            PartitionClockModel pcm = new PartitionClockModel(this, DEFAULT_NAME, partition, partition.getPartitionTreeModel());
            partition.setPartitionClockModel(pcm);
        }
        // if (partition.getPartitionSubstitutionModel() == null) {
        //     PartitionSubstitutionModel substModel = new PartitionSubstitutionModel(this, partition.getName(),
        //             partition);
        //     partition.setPartitionSubstitutionModel(substModel);
        // }

        ContinuousComponentOptions comp = (ContinuousComponentOptions) getComponentOptions(ContinuousComponentOptions.class);
        comp.createParameters(this);

        return selRow;
    }

    // ++++++++++++++++++++ message bar +++++++++++++++++

    public String statusMessage() {
        // String message = "<html><p>";
        String message = "";
        if (hasData()) {
            message += "Data: " + taxonList.getTaxonCount() + " taxa, ";
            message += dataPartitions.size() + (dataPartitions.size() > 1 ? " partitions" : " partition");

            if (userTrees.size() > 0) {
                message += ", " + userTrees.size() + " user" +
                        (userTrees.size() > 1 ? " trees" : " tree");
            }

            // if (hasPhylogeographic()) {
            //      message += ";    Phylogeographic Analysis";
            // }

            // message += "; " + clockModelOptions.statusMessageClockModel();

        } else if (userTrees.size() > 0) { // TODO
            message += "Trees only : " + userTrees.size() +
                    (userTrees.size() > 1 ? " trees, " : " tree, ") +
                    taxonList.getTaxonCount() + " taxa";
        } else if (taxonList != null && taxonList.getTaxonCount() > 0) {  // TODO
            message += "Taxa only: " + taxonList.getTaxonCount() + " taxa";
        } else {
            message += "No data loaded - select 'Import Data...' from the 'File' menu.";
        }
        // message += "</p></html>";
        return message;
    }

    public List<Object> getKeysFromValue(Map<?, ?> hm, Object value) {
        List<Object> list = new ArrayList<Object>();
        for (Object o : hm.keySet()) {
            if (hm.get(o).equals(value)) {
                list.add(o);
            }
        }
        return list;
    }

    public Taxa getTaxa(String taxaName) {
        for (Taxa taxa : taxonSets) {
            if (taxa.getId().equalsIgnoreCase(taxaName)) return taxa;
        }
        return null;
    }

    public int getTaxaIndex(String taxaName) {
        for (int i = 0; i < taxonSets.size(); i++) {
            if (taxonSets.get(i).getId().equalsIgnoreCase(taxaName)) return i;
        }
        return -1;
    }

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Data options
    public Taxa taxonList = null; // union set of all taxa in all partitions. todo change to List<Taxa> regarding data type?

    public List<Taxa> taxonSets = new ArrayList<Taxa>();
    public Map<Taxa, Boolean> taxonSetsMono = new HashMap<Taxa, Boolean>();
    public Map<Taxa, Double> taxonSetsHeights = new HashMap<Taxa, Double>();
    public Map<Taxa, Boolean> taxonSetsIncludeStem = new HashMap<Taxa, Boolean>();
    public Map<Taxa, PartitionTreeModel> taxonSetsTreeModel = new HashMap<Taxa, PartitionTreeModel>();

    public boolean useTipDates = false;

    public DateUnitsType datesUnits = DateUnitsType.YEARS;
    public DateUnitsType datesDirection = DateUnitsType.FORWARDS;
    public double maximumTipHeight = 0.0;
    public int translation = 0;

    public Date originDate = null;


    public DateGuesser dateGuesser = new DateGuesser();
    // public TraitGuesser traitGuesser = new TraitGuesser();
    //
    // public List<String> selecetedTraits = new ArrayList<String>();
    // public Map<String, TraitGuesser.TraitType> traitTypes = new HashMap<String, TraitGuesser.TraitType>();

    // Data
    public List<AbstractPartitionData> dataPartitions = new ArrayList<AbstractPartitionData>();
    public List<TraitData> traits = new ArrayList<TraitData>();

    public List<List<PartitionData>> multiPartitionLists = new ArrayList<List<PartitionData>>();
    public List<AbstractPartitionData> otherPartitions = new ArrayList<AbstractPartitionData>();

    // list of starting tree from user import
    public Map<String, TreeHolder> userTrees = new HashMap<>();

    public boolean unlinkPartitionRates = true;

    public Units.Type units = Units.Type.YEARS;

    // Operator schedule options
    public OperatorSchedule.OptimizationTransform optimizationTransform = OperatorSchedule.DEFAULT_TRANSFORM;

    // MCMC options
    public int chainLength = 10000000;
    public int logEvery = 1000;
    public int echoEvery = 1000;
    public int checkpointEvery = 1000000;
    public int burnIn = 100000;
    public String fileName = null;
    public boolean autoOptimize = true;
    public boolean performTraceAnalysis = false;
    public boolean generateCSV = true;  // until/if a button
    public boolean samplePriorOnly = false;

    public String fileNameStem = MCMCPanel.DEFAULT_FILE_NAME_STEM;
    public String logFileName = null;
    public String checkpointFileName = null;

    public boolean generateDemographicLogFile = false;
    public String demographicModelName = null;
    public String demographicLogFileName = null;

    public boolean allowOverwriteLog = false;
    // public boolean mapTreeLog = false;
    // public String mapTreeFileName = null;
    public List<String> treeFileName = new ArrayList<String>();
    public boolean substTreeLog = false;
    public List<String> substTreeFileName = new ArrayList<String>();
    public boolean operatorAnalysis = true;
    public String operatorAnalysisFileName = null;

    public GlobalModelOptions globalModelOptions = new GlobalModelOptions(this);
    public ClockModelOptions clockModelOptions = new ClockModelOptions(this);
    public TreeModelOptions treeModelOptions = new TreeModelOptions(this);

    public boolean useClassicOperatorsAndPriors = false;

    public OperatorSetType operatorSetType = OperatorSetType.DEFAULT;

    public boolean useClassicOperatorsAndPriors() {
        return useClassicOperatorsAndPriors;
    }

    public boolean useNuRelativeRates() {
        return !useClassicOperatorsAndPriors() || !NEW_RELATIVE_RATE_PARAMETERIZATION;
    }

    public boolean useGammaPriorPopSize() {
        return !useClassicOperatorsAndPriors() || !POPSIZE_GAMMA_PRIOR;
    }

    public boolean useNewGTR() {
        return !useClassicOperatorsAndPriors() || !NEW_GTR_PARAMETERIZATION;
    }

    public boolean usePInvRandomWalk() {
        return !useClassicOperatorsAndPriors() || !LOGIT_PINV_KERNEL;
    }

    public boolean useNewFrequenciesPrior() {
        return !useClassicOperatorsAndPriors() || !FREQUENCIES_DIRICHLET_PRIOR;
    }

    public boolean logCoalescentEventsStatistic = false;

}