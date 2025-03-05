/*
 * BeastCheckpointer.java
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

package dr.app.checkpoint;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.markovchain.MarkovChain;
import dr.inference.markovchain.MarkovChainListener;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorSchedule;
import dr.inference.state.*;
import dr.math.MathUtils;

import java.io.*;
import java.text.SimpleDateFormat;
import java.time.LocalTime;
import java.time.format.DateTimeFormatter;
import java.util.*;

/**
 * A state loader / saver
 * @author Andrew Rambaut
 * @author Guy Baele
 */
public class BeastCheckpointer implements StateLoaderSaver {

    private static BeastCheckpointer single_instance = null;

    private static final boolean DEBUG = false;

    // A debugging flag to do a check that the state gives the same likelihood after loading
    private static final boolean CHECK_LOAD_STATE = true;

    public final static String LOAD_STATE_FILE = "load.state.file";
    public final static String SAVE_STATE_FILE = "save.state.file";
    public final static String SAVE_STATE_AT = "save.state.at";
    public final static String SAVE_STATE_EVERY = "save.state.every";
    public final static String SAVE_STEM = "save.state.stem";
    public final static String SAVE_STATE_TIME = "save.state.time";

    public final static String FORCE_RESUME = "force.resume";
    public final static String CHECKPOINT_SEED = "checkpoint.seed";

    public final static String FULL_CHECKPOINT_PRECISION = "full.checkpoint.precision";

    private String loadStateFileName;
    private String saveStateFileName;
    private String stemFileName;

    private boolean forceResume = false;

    private boolean useFullPrecision;

    private final List<MarkovChainListener> listeners = new ArrayList<MarkovChainListener>();

    public static synchronized BeastCheckpointer getInstance(String checkpointFileName, int checkpointEvery, int checkpointFinal, boolean overwrite) {
        if (single_instance == null) {
            single_instance = new BeastCheckpointer(checkpointFileName, checkpointEvery, checkpointFinal, overwrite);
        }

        return single_instance;
    }

    private BeastCheckpointer(String checkpointFileName, int checkpointEvery, int checkpointFinal, boolean overwrite) {

        if (checkpointFileName == null && checkpointEvery == -1 && checkpointFinal == -1) {

            loadStateFileName = System.getProperty(LOAD_STATE_FILE, null);
            saveStateFileName = System.getProperty(SAVE_STATE_FILE, null);

            stemFileName = System.getProperty(SAVE_STEM, null);

            if (System.getProperty(SAVE_STATE_AT) != null) {
                final long saveStateAt = Long.parseLong(System.getProperty(SAVE_STATE_AT));
                listeners.add(new StateSaverChainListener(BeastCheckpointer.this, saveStateAt,false));
            }
            if (System.getProperty(SAVE_STATE_EVERY) != null) {
                final long saveStateEvery = Long.parseLong(System.getProperty(SAVE_STATE_EVERY));
                listeners.add(new StateSaverChainListener(BeastCheckpointer.this, saveStateEvery,true));
            }

        } else {

            loadStateFileName = System.getProperty(LOAD_STATE_FILE, null);
            saveStateFileName = checkpointFileName;

            stemFileName = System.getProperty(SAVE_STEM, null);

            listeners.add(new StateSaverChainListener(BeastCheckpointer.this, checkpointFinal,false));
            listeners.add(new StateSaverChainListener(BeastCheckpointer.this, checkpointEvery,true));

        }
        if (System.getProperty(SAVE_STATE_TIME) != null) {
            LocalTime saveTime = LocalTime.parse(System.getProperty(SAVE_STATE_TIME),
                    DateTimeFormatter.ofPattern("HH:mm:ss"));
            int saveSeconds = saveTime.toSecondOfDay();
            listeners.add(new TimedStateSaverChainListener(BeastCheckpointer.this, saveSeconds));
        }

        useFullPrecision = (System.getProperty(FULL_CHECKPOINT_PRECISION) != null) &&
                System.getProperty(FULL_CHECKPOINT_PRECISION).equalsIgnoreCase("true");

        Factory.INSTANCE = new Factory() {
            @Override
            public StateLoader getInitialStateLoader() {
                if (loadStateFileName == null) {
                    return null;
                } else {
                    return getStateLoaderObject();
                }
            }

            @Override
            public MarkovChainListener[] getStateSaverChainListeners() {
                return listeners.toArray(new MarkovChainListener[0]);
            }

            @Override
            public StateLoaderSaver getStateLoaderSaver(final File loadFile, final File saveFile) {
                return new StateLoaderSaver() {

                    @Override
                    public boolean saveState(MarkovChain markovChain, long state, double lnL) {
                        return BeastCheckpointer.this.writeStateToFile(saveFile, state, lnL, markovChain);
                    }

                    @Override
                    public long loadState(MarkovChain markovChain, double[] savedLnL) {
                        return BeastCheckpointer.this.readStateFromFile(loadFile, markovChain, savedLnL);
                    }

                    @Override
                    public void checkLoadState(double savedLnL, double lnL) {
                        // do nothing.
                    }
                };
            }
        };

    }

    //only here for reasons of inheritance (of the CheckPointModifier class)
    protected BeastCheckpointer() {

    }

    private BeastCheckpointer getStateLoaderObject() {
        return this;
    }

    @Override
    public boolean saveState(MarkovChain markovChain, long state, double lnL) {
        String fileName = "";
        if (stemFileName != null && this.saveStateFileName == null) {
            fileName = stemFileName + "_" + state;
        } else if (stemFileName != null) {
            fileName = stemFileName + "_" + this.saveStateFileName;
        } else {
            String timeStamp = new SimpleDateFormat("yyyy.MM.dd.HH.mm.ss").format(Calendar.getInstance().getTime());
            fileName = (this.saveStateFileName != null ? this.saveStateFileName : "beast_state_" + timeStamp);
        }
        return writeStateToFile(new File(fileName), state, lnL, markovChain);
    }

    @Override
    public long loadState(MarkovChain markovChain, double[] savedLnL) {
        return readStateFromFile(new File(loadStateFileName), markovChain, savedLnL);
    }

    @Override
    public void checkLoadState(double savedLnL, double lnL) {

        if (CHECK_LOAD_STATE) {

            if (System.getProperty(FORCE_RESUME) != null) {
                forceResume = Boolean.parseBoolean(System.getProperty(FORCE_RESUME));
                //System.out.println("FORCE_RESUME? " + System.getProperty(FORCE_RESUME) + " (" + forceResume + "");
            }

            //first perform a simple check for equality of two doubles
            //when this test fails, go over the digits
            if (forceResume) {
                System.out.println("Forcing analysis to resume regardless of recomputed likelihood values ("  + lnL + " vs. " + savedLnL + ").");
            } else if (lnL != savedLnL) {

                System.out.println("COMPARING LIKELIHOODS: " + lnL + " vs. " + savedLnL);

                //15 is the floor value for the number of decimal digits when representing a double
                //checking for 15 identical digits below
                String originalString = Double.toString(savedLnL);
                String restoredString = Double.toString(lnL);
                System.out.println(lnL + "    " + originalString);
                System.out.println(savedLnL + "    " + restoredString);
                //assume values will be nearly identical
                int digits = 0;
                for (int i = 0; i < Math.min(originalString.length(), restoredString.length()); i++) {
                    if (originalString.charAt(i) == restoredString.charAt(i)) {
                        if (!(originalString.charAt(i) == '-' || originalString.charAt(i) == '.')) {
                            digits++;
                        }
                    } else {
                        break;
                    }
                }
                //System.err.println("digits = " + digits);

                //the translation table is constructed upon loading the initial tree for the analysis
                //if that tree is user-defined or a UPGMA starting tree, the starting seed won't matter
                //if that tree is constructed at random, we currently don't provide a way to retrieve the original translation table

                if (digits < 15) {

                    //currently use the general BEAST -threshold argument
                    //TODO Evaluate whether a checkpoint-specific threshold option is required or useful
                    double threshold = 1E-10;
                    if (System.getProperty("mcmc.evaluation.threshold") != null) {
                        threshold = Double.parseDouble(System.getProperty("mcmc.evaluation.threshold"));
                    }
                    if (Math.abs(lnL - savedLnL) > threshold) {
                        throw new RuntimeException("Saved lnL does not match recomputed value for loaded state: stored lnL: " + savedLnL +
                                ", recomputed lnL: " + lnL + " (difference " + (savedLnL - lnL) + ")." +
                                "\nYour XML may require the construction of a randomly generated starting tree. " +
                                "Try resuming the analysis by using the same starting seed as for the original BEAST run.");
                    } else {
                        System.out.println("Saved lnL does not match recomputed value for loaded state: stored lnL: " + savedLnL +
                                ", recomputed lnL: " + lnL + " (difference " + (savedLnL - lnL) + ")." +
                                "\nThreshold of " + threshold + " for restarting analysis not exceeded; continuing ...");
                    }
                }

            } else {
                System.out.println("IDENTICAL LIKELIHOODS");
                System.out.println("lnL" + " = " + lnL);
                System.out.println("savedLnL[0]" + " = " + savedLnL);
            }

        }
    }

    protected boolean writeStateToFile(File file, long state, double lnL, MarkovChain markovChain) {
        OperatorSchedule operatorSchedule = markovChain.getSchedule();

        OutputStream fileOut = null;
        try {
            fileOut = new FileOutputStream(file);
            PrintStream out = useFullPrecision ? new CheckpointPrintStream(fileOut) : new PrintStream(fileOut);

            ArrayList<TreeParameterModel> traitModels = new ArrayList<TreeParameterModel>();

            int[] rngState = MathUtils.getRandomState();
            out.print("rng");
            for (int i = 0; i < rngState.length; i++) {
                out.print("\t");
                out.print(rngState[i]);
            }
            out.println();

            out.print("state\t");
            out.println(state);

            out.print("lnL\t");
            out.println(lnL);

            for (Parameter parameter : Parameter.CONNECTED_PARAMETER_SET) {
                if (!parameter.isImmutable()) {
                    out.print("parameter");
                    out.print("\t");
                    out.print(parameter.getParameterName());
                    out.print("\t");
                    out.print(parameter.getDimension());
                    for (int dim = 0; dim < parameter.getDimension(); dim++) {
                        out.print("\t");
                        out.print(parameter.getParameterUntransformedValue(dim));
                    }
                    out.print("\n");
                }
            }

            for (int i = 0; i < operatorSchedule.getOperatorCount(); i++) {
                MCMCOperator operator = operatorSchedule.getOperator(i);
                out.print("operator");
                out.print("\t");
                out.print(operator.getOperatorName());
                out.print("\t");
                out.print(operator.getAcceptCount());
                out.print("\t");
                out.print(operator.getRejectCount());
                if (operator instanceof AdaptableMCMCOperator) {
                    out.print("\t");
                    out.print(((AdaptableMCMCOperator)operator).getAdaptableParameter());
                    out.print("\t");
                    out.print(((AdaptableMCMCOperator)operator).getAdaptationCount());
                }
                out.println();
            }

            //check up front if there are any TreeParameterModel objects
            for (Model model : Model.CONNECTED_MODEL_SET) {
                if (model instanceof TreeParameterModel) {
                    if (DEBUG) {
                        System.out.println("\nSave TreeParameterModel: " + model.getClass().getSimpleName());
                    }
                    traitModels.add((TreeParameterModel) model);
                }
            }

            for (Model model : Model.CONNECTED_MODEL_SET) {

                if (model instanceof TreeModel) {
                    out.print("tree");
                    out.print("\t");
                    out.println(model.getModelName());

                    //replace Newick format by printing general graph structure
                    //out.println(((TreeModel) model).getNewick());

                    out.println("#node height taxon");
                    int nodeCount = ((TreeModel) model).getNodeCount();
                    out.println(nodeCount);
                    for (int i = 0; i < nodeCount; i++) {
                        out.print(((TreeModel) model).getNode(i).getNumber());
                        out.print("\t");
                        out.print(((TreeModel) model).getNodeHeight(((TreeModel) model).getNode(i)));
                        if (((TreeModel) model).isExternal(((TreeModel) model).getNode(i))) {
                            out.print("\t");
                            out.print(((TreeModel) model).getNodeTaxon(((TreeModel) model).getNode(i)).getId());
                        }
                        out.println();
                    }

                    out.println("#edges");
                    out.println("#child-node parent-node L/R-child traits");

                    out.println(nodeCount);
                    for (int i = 0; i < nodeCount; i++) {
                        NodeRef parent = ((TreeModel) model).getParent(((TreeModel) model).getNode(i));
                        if (parent != null) {
                            out.print(((TreeModel) model).getNode(i).getNumber());
                            out.print("\t");
                            out.print(((TreeModel) model).getParent(((TreeModel) model).getNode(i)).getNumber());
                            out.print("\t");

                            if ((((TreeModel) model).getChild(parent, 0) == ((TreeModel) model).getNode(i))) {
                                //left child
                                out.print(0);
                            } else if ((((TreeModel) model).getChild(parent, 1) == ((TreeModel) model).getNode(i))) {
                                //right child
                                out.print(1);
                            } else {
                                throw new RuntimeException("Operation currently only supported for nodes with 2 children.");
                            }

                            //only print the TreeParameterModel that matches the TreeModel currently being written
                            for (TreeParameterModel tpm : traitModels) {
                                if (model == tpm.getTreeModel()) {
                                    out.print("\t");
                                    out.print(tpm.getNodeValue((TreeModel) model, ((TreeModel) model).getNode(i)));
                                }
                            }
                            out.println();
                        } else {
                            if (DEBUG) {
                                System.out.println(((TreeModel) model).getNode(i) + " has no parent.");
                            }
                        }
                    }

                }

            }

            out.close();
            fileOut.close();
        } catch (IOException ioe) {
            System.err.println("Unable to write file: " + ioe.getMessage());
            return false;
        }

        if (DEBUG) {
            for (Likelihood likelihood : Likelihood.CONNECTED_LIKELIHOOD_SET) {
                System.err.println(likelihood.getId() + ": " + likelihood.getLogLikelihood());
            }
        }

        return true;
    }

    protected long readStateFromFile(File file, MarkovChain markovChain, double[] lnL) {

        DoubleParser parser = useFullPrecision ? DoubleParser.HEX : DoubleParser.TEXT;

        OperatorSchedule operatorSchedule = markovChain.getSchedule();

        long state = -1;

        ArrayList<TreeParameterModel> traitModels = new ArrayList<TreeParameterModel>();

        try {
            FileReader fileIn = new FileReader(file);
            BufferedReader in = new BufferedReader(fileIn);

            int[] rngState = null;

            String line = in.readLine();
            String[] fields = line.split("\t");
            if (fields[0].equals("rng")) {
                // if there is a random number generator state present then load it...
                try {
                    rngState = new int[fields.length - 1];
                    for (int i = 0; i < rngState.length; i++) {
                        rngState[i] = Integer.parseInt(fields[i + 1]);
                    }

                } catch (NumberFormatException nfe) {
                    throw new RuntimeException("Unable to read state number from state file");
                }

                line = in.readLine();
                fields = line.split("\t");
            }

            try {
                if (!fields[0].equals("state")) {
                    throw new RuntimeException("Unable to read state number from state file");
                }
                state = Long.parseLong(fields[1]);
            } catch (NumberFormatException nfe) {
                throw new RuntimeException("Unable to read state number from state file");
            }

            line = in.readLine();
            fields = line.split("\t");
            try {
                if (!fields[0].equals("lnL")) {
                    throw new RuntimeException("Unable to read lnL from state file");
                }
                if (lnL != null) {
                    lnL[0] = parser.parseDouble(fields[1]);
                }
            } catch (NumberFormatException nfe) {
                throw new RuntimeException("Unable to read lnL from state file");
            }

            for (Parameter parameter : Parameter.CONNECTED_PARAMETER_SET) {

                if (!parameter.isImmutable()) {
                    line = in.readLine();
                    fields = line.split("\t");
                    //if (!fields[0].equals(parameter.getParameterName())) {
                    //  System.err.println("Unable to match state parameter: " + fields[0] + ", expecting " + parameter.getParameterName());
                    //}
                    int dimension = Integer.parseInt(fields[2]);

                    if (dimension != parameter.getDimension()) {
                        System.err.println("Unable to match state parameter dimension: " + dimension + ", expecting " + parameter.getDimension() + " for parameter: " + parameter.getParameterName());
                        System.err.print("Read from file: ");
                        for (int i = 0; i < fields.length; i++) {
                            System.err.print(fields[i] + "\t");
                        }
                        System.err.println();
                    }

                    if (fields[1].equals("branchRates.categories.rootNodeNumber")) {
                        // System.out.println("eek");
                        double value = parser.parseDouble(fields[3]);
                        parameter.setParameterValue(0, value);
                        if (DEBUG) {
                            System.out.println("restoring " + fields[1] + " with value " + value);
                        }
                    } else {
                        if (DEBUG) {
                            System.out.print("restoring " + fields[1] + " with values ");
                        }
                        for (int dim = 0; dim < parameter.getDimension(); dim++) {
                            try {
                                parameter.setParameterUntransformedValue(dim, parser.parseDouble(fields[dim + 3]));
                            } catch (RuntimeException rte) {
                                System.err.println(rte);
                                continue;
                            }
                            if (DEBUG) {
                                System.out.print(parser.parseDouble(fields[dim + 3]) + " ");
                            }
                        }
                        if (DEBUG) {
                            System.out.println();
                        }
                    }
                }
            }

            for (int i = 0; i < operatorSchedule.getOperatorCount(); i++) {
                //TODO we can no longer assume these are in the right order
                //TODO best parse all the "operator" lines and store them so we can mix and match within this for loop
                //TODO does not only apply to the operators but also to the parameters
                //TODO test using additional tip-date sampling compared to previous run
                MCMCOperator operator = operatorSchedule.getOperator(i);
                line = in.readLine();
                fields = line.split("\t");
                if (!fields[1].equals(operator.getOperatorName())) {
                    throw new RuntimeException("Unable to match " + operator.getOperatorName() + " operator: " + fields[1]);
                }
                if (fields.length < 4) {
                    throw new RuntimeException("Operator missing values: " + fields[1]);
                }
                operator.setAcceptCount(Integer.parseInt(fields[2]));
                operator.setRejectCount(Integer.parseInt(fields[3]));
                if (operator instanceof AdaptableMCMCOperator) {
                    if (fields.length != 6) {
                        throw new RuntimeException("Coercable operator missing parameter: " + fields[1]);
                    }
                    ((AdaptableMCMCOperator)operator).setAdaptableParameter(parser.parseDouble(fields[4]));
                    ((AdaptableMCMCOperator)operator).setAdaptationCount(Long.parseLong(fields[5]));
                }
            }

            // load the tree models last as we get the node heights from the tree (not the parameters which
            // which may not be associated with the right node
            Set<String> expectedTreeModelNames = new LinkedHashSet<>();

            //store list of TreeModels for debugging purposes
            ArrayList<TreeModel> treeModelList = new ArrayList<TreeModel>();

            for (Model model : Model.CONNECTED_MODEL_SET) {

                if (model instanceof TreeModel) {
                    if (DEBUG) {
                        System.out.println("model " + model.getModelName());
                    }
                    treeModelList.add((TreeModel)model);
                    expectedTreeModelNames.add(model.getModelName());
                    if (DEBUG) {
                        System.out.println("\nexpectedTreeModelNames:");
                        for (String s : expectedTreeModelNames) {
                            System.out.println(s);
                        }
                        System.out.println();
                    }
                } else {
                    if (DEBUG) {
                        System.out.println("Not a TreeModel: " + model.getModelName());
                    }
                }

                //first add all TreeParameterModels to a list
                if (model instanceof TreeParameterModel) {
                    if (DEBUG) {
                        System.out.println("\nLoad TreeParameterModel: " + model.getClass().getSimpleName());
                    }
                    traitModels.add((TreeParameterModel)model);
                    if (DEBUG) {
                        System.out.println("TreeParameterModel: " + model.getModelName());
                    }
                }

            }

            //explicitly link TreeModel (using its unique ID) to a list of TreeParameterModels
            //this information is currently not yet used
            HashMap<String, ArrayList<TreeParameterModel>> linkedModels = new HashMap<String, ArrayList<TreeParameterModel>>();
            for (String name : expectedTreeModelNames) {
                ArrayList<TreeParameterModel> tpmList = new ArrayList<TreeParameterModel>();
                for (TreeParameterModel tpm : traitModels) {
                    if (tpm.getTreeModel().getId().equals(name)) {
                        tpmList.add(tpm);
                        if (DEBUG) {
                            System.out.println("TreeModel: " + name + " has been assigned TreeParameterModel: " + tpm.toString());
                        }
                    }
                }
                linkedModels.put(name, tpmList);
            }

            line = in.readLine();
            fields = line.split("\t");
            // Read in all (possibly more than one) trees
            while (fields[0].equals("tree")) {

                if (DEBUG) {
                    System.out.println("\ntree: " + fields[1]);
                }

                for (Model model : Model.CONNECTED_MODEL_SET) {
                    if (model instanceof TreeModel && fields[1].equals(model.getModelName())) {
                        line = in.readLine();
                        line = in.readLine();
                        fields = line.split("\t");
                        //read number of nodes
                        int nodeCount = Integer.parseInt(fields[0]);
                        double[] nodeHeights = new double[nodeCount];
                        String[] taxaNames = new String[(nodeCount+1)/2];

                        for (int i = 0; i < nodeCount; i++) {
                            line = in.readLine();
                            fields = line.split("\t");
                            nodeHeights[i] = parser.parseDouble(fields[1]);
                            if (i < taxaNames.length) {
                                taxaNames[i] = fields[2];
                            }
                        }

                        //on to reading edge information
                        line = in.readLine();
                        line = in.readLine();
                        line = in.readLine();
                        fields = line.split("\t");

                        int edgeCount = Integer.parseInt(fields[0]);
                        if (DEBUG) {
                            System.out.println("edge count = " + edgeCount);
                            System.out.println("model: " + model.getId());
                            System.out.println("linkedModels size = " + linkedModels.size());
                            System.out.println(linkedModels.get(model.getId()));
                        }

                        //create data matrix of doubles to store information from list of TreeParameterModels
                        //size of matrix depends on the number of TreeParameterModels assigned to a TreeModel
                        double[][] traitValues = new double[linkedModels.get(model.getId()).size()][edgeCount];

                        //create array to store whether a node is left or right child of its parent
                        //can be important for certain tree transition kernels
                        int[] childOrder = new int[edgeCount];
                        for (int i = 0; i < childOrder.length; i++) {
                            childOrder[i] = -1;
                        }

                        int[] parents = new int[edgeCount];
                        for (int i = 0; i < edgeCount; i++){
                            parents[i] = -1;
                        }
                        for (int i = 0; i < edgeCount-1; i++) {
                            line = in.readLine();
                            if (line != null) {
                                if (DEBUG) {
                                    System.out.println("DEBUG: " + line);
                                }
                                fields = line.split("\t");
                                parents[Integer.parseInt(fields[0])] = Integer.parseInt(fields[1]);
                                // childOrder[i] = Integer.parseInt(fields[2]);
                                childOrder[Integer.parseInt(fields[0])] = Integer.parseInt(fields[2]);
                                for (int j = 0; j < linkedModels.get(model.getId()).size(); j++) {
                                    //   traitValues[j][i] = parser.parseDouble(fields[3+j]);
                                    traitValues[j][Integer.parseInt(fields[0])] = parser.parseDouble(fields[3+j]);
                                }
                            }
                        }

                        //perform magic with the acquired information
                        if (DEBUG) {
                            System.out.println("adopting tree structure");
                        }

                        //adopt the loaded tree structure
                        ((TreeModel) model).beginTreeEdit();
                        ((TreeModel) model).adoptTreeStructure(parents, nodeHeights, childOrder, taxaNames);
                        if (traitModels.size() > 0) {
                            System.out.println("adopting " + traitModels.size() + " trait models to treeModel " + ((TreeModel)model).getId());
                            ((TreeModel) model).adoptTraitData(parents, traitModels, traitValues, taxaNames);
                        }
                        ((TreeModel) model).endTreeEdit();

                        expectedTreeModelNames.remove(model.getModelName());

                    }

                }

                line = in.readLine();
                if (line != null) {
                    fields = line.split("\t");
                }

            }

            if (expectedTreeModelNames.size() > 0) {
                StringBuilder sb = new StringBuilder();
                for (String notFoundName : expectedTreeModelNames) {
                    sb.append("Expecting, but unable to match state parameter:" + notFoundName + "\n");
                }
                throw new RuntimeException("\n" + sb.toString());
            }

            if (DEBUG) {
                System.out.println("\nDouble checking:");
                for (Parameter parameter : Parameter.CONNECTED_PARAMETER_SET) {
                    if (parameter.getParameterName() != null && parameter.getParameterName().equals("branchRates.categories.rootNodeNumber")) {
                        System.out.println(parameter.getParameterName() + ": " + parameter.getParameterValue(0));
                    }
                }
                System.out.println("\nPrinting trees:");
                for (TreeModel tm : treeModelList) {
                    System.out.println(tm.getId() + ": ");
                    System.out.println(tm.getNewick());
                }
            }

            if (System.getProperty(BeastCheckpointer.CHECKPOINT_SEED) != null) {
                MathUtils.setSeed(Long.parseLong(System.getProperty(BeastCheckpointer.CHECKPOINT_SEED)));
            } else if (rngState != null) {
                MathUtils.setRandomState(rngState);
            }

            in.close();
            fileIn.close();

            //This shouldn't be necessary and if it is then it might be hiding a bug...
            /*for (Likelihood likelihood : Likelihood.CONNECTED_LIKELIHOOD_SET) {
                likelihood.makeDirty();
            }*/

        } catch (IOException ioe) {
            throw new RuntimeException("Unable to read file: " + ioe.getMessage());
        }

        return state;
    }

    class CheckpointPrintStream extends PrintStream {
        public CheckpointPrintStream(OutputStream out) {
            super(out);
        }

        public void print(double d) {
            String byteString = Long.toHexString(Double.doubleToRawLongBits(d));
            super.print(d + "/" + byteString);
        }

        public void println(double d) {
            String byteString = Long.toHexString(Double.doubleToRawLongBits(d));
            super.println(d + "/" + byteString);
        }
    }

    enum DoubleParser {
        TEXT {
            public double parseDouble(String string) {
                return Double.parseDouble(string);
            }
        },
        HEX {
            public double parseDouble(String string) {
                String[] strings = string.split("/");
                return convertHexStrToDouble(strings[1]);
            }

            private double convertHexStrToDouble(String input) {
                // convert the input to positive, as needed
                String s2 = preprocess(input);
                boolean negative = true;
                // if the original equals the new string, then it is not negative
                if (input.equalsIgnoreCase(s2)) {
                    negative = false;
                }

                // convert the hex string to long
                long doubleAsLongReverse = Long.parseLong(s2, 16);

                // Convert the long back into the original double
                double doubleOutput = Double.longBitsToDouble(doubleAsLongReverse);

                // return as a negative value, as needed
                if (negative) {
                    return -doubleOutput;
                } else {
                    return doubleOutput;
                }
            }

            private String preprocess(String doubleAsHexString) {
                // get the first char and convert it to an int
                String s0 = doubleAsHexString.substring(0, 1);
                int int1 = Integer.parseInt(s0, 16);

                // if the int is < 8, then the string is not negative
                // and is returned without further processing
                if (int1 < 8) {
                    return doubleAsHexString;
                }

                // otherwise subtract 8
                int1 = int1 - 8;
                s0 = Integer.toString(int1);

                // don't prepend a "0"
                if (int1 == 0) {
                    s0 = "";
                }

                // return the string with a new inital char
                return s0 + doubleAsHexString.substring(1);
            }

        };

        public abstract double parseDouble(String string);
    }

}
