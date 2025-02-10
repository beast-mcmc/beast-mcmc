/*
 * MarkovChain.java
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

package dr.inference.markovchain;

import dr.evomodel.continuous.GibbsIndependentCoalescentOperator;
import dr.inference.model.*;
import dr.inference.operators.*;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 * A concrete markov chain. This is final as the only things that should need
 * overriding are in the delegates (prior, likelihood, schedule and acceptor).
 * The design of this class is to be fairly immutable as far as settings goes.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public final class MarkovChain implements Serializable {
    private static final long serialVersionUID = 181L;

    private final static boolean DEBUG = false;
    private final static boolean PROFILE = true;

    public static final double EVALUATION_TEST_THRESHOLD = 1e-1;

    private final OperatorSchedule schedule;
    private final Acceptor acceptor;
    private final Likelihood likelihood;

    private boolean pleaseStop = false;
    private boolean isStopped = false;
    private double bestScore, currentScore, initialScore;
    private long currentLength;

    private final boolean useAdaptation;
    private boolean isCurrentlyAdapting;

    private final boolean useSmoothedAcceptanceProbability;


    private final long fullEvaluationCount;
    private final int minOperatorCountForFullEvaluation;

    private double evaluationTestThreshold = EVALUATION_TEST_THRESHOLD;


    public MarkovChain(Likelihood likelihood,
                       OperatorSchedule schedule, Acceptor acceptor,
                       long fullEvaluationCount, int minOperatorCountForFullEvaluation, double evaluationTestThreshold,
                       boolean useAdaptation, boolean useSmoothedAcceptanceProbability) {

        currentLength = 0;
        this.likelihood = likelihood;
        this.schedule = schedule;
        this.acceptor = acceptor;
        this.useAdaptation = useAdaptation;
        this.isCurrentlyAdapting = useAdaptation;
        this.useSmoothedAcceptanceProbability = useSmoothedAcceptanceProbability;

        this.fullEvaluationCount = fullEvaluationCount;
        this.minOperatorCountForFullEvaluation = minOperatorCountForFullEvaluation;
        this.evaluationTestThreshold = evaluationTestThreshold;

        Likelihood.CONNECTED_LIKELIHOOD_SET.add(likelihood);
        Likelihood.CONNECTED_LIKELIHOOD_SET.addAll(likelihood.getLikelihoodSet());

        for (Likelihood l : Likelihood.FULL_LIKELIHOOD_SET) {
            if (!Likelihood.CONNECTED_LIKELIHOOD_SET.contains(l)) {
                System.err.println("WARNING: Likelihood component, " + l.getId() + ", created but not used in the MCMC");
            }
        }

        currentScore = evaluate(likelihood);
    }

    /**
     * Resets the markov chain
     */
    public void reset() {
        currentLength = 0;

        // reset operator acceptance levels
        for (int i = 0; i < schedule.getOperatorCount(); i++) {
            schedule.getOperator(i).reset();
        }
    }

    /**
     * Run the chain for a given number of states.
     *
     * @param length number of states to run the chain.
     */
    public long runChain(long length, boolean disableAdaptation) {

        isCurrentlyAdapting = useAdaptation && !disableAdaptation;

        likelihood.makeDirty();
        currentScore = evaluate(likelihood);

        long currentState = currentLength;

        final Model currentModel = likelihood.getModel();

        if (currentState == 0) {
            initialScore = currentScore;
            bestScore = currentScore;
            fireBestModel(currentState, currentModel);
        }

        if (currentScore == Double.NEGATIVE_INFINITY) {

            // identify which component of the score is zero...
            String message = "The initial likelihood is zero";
            if (likelihood instanceof CompoundLikelihood) {
                message += ": " + ((CompoundLikelihood) likelihood).getDiagnosis();
            } else if (likelihood instanceof PathLikelihood) {
                message += ": " + ((CompoundLikelihood)((PathLikelihood) likelihood).getSourceLikelihood()).getDiagnosis();
                message += ": " + ((CompoundLikelihood)((PathLikelihood) likelihood).getDestinationLikelihood()).getDiagnosis();
            } else {
                message += ".";
            }
            throw new IllegalArgumentException(message);
        } else if (currentScore == Double.POSITIVE_INFINITY || Double.isNaN(currentScore)) {
            String message = "A likelihood returned with a numerical error (" + currentScore + ")";
            if (likelihood instanceof CompoundLikelihood) {
                message += ": " + ((CompoundLikelihood) likelihood).getDiagnosis();
            } else {
                message += ".";
                Set<Likelihood> lSet = likelihood.getLikelihoodSet();
                for (Likelihood l : lSet) {
                    message += "\n  " + l.prettyName() + " = " + l.getLogLikelihood();
                }
            }
            throw new IllegalArgumentException(message);
        }

        pleaseStop = false;
        isStopped = false;

        //int otfcounter = onTheFlyOperatorWeights > 0 ? onTheFlyOperatorWeights : 0;

        double[] logr = {0.0};

        boolean usingFullEvaluation = true;
        // set ops count in mcmc element instead
        if (fullEvaluationCount == 0) // Temporary solution until full code review
            usingFullEvaluation = false;
        boolean fullEvaluationError = false;

        Map<String, Double> diagnosticDensities = null;

        while (!pleaseStop && (currentState < (currentLength + length))) {

            // periodically log states
            fireCurrentModel(currentState, currentModel);

            if (pleaseStop) {
                isStopped = true;
                break;
            }

            // Get the operator
            final int op = schedule.getNextOperatorIndex();
            final MCMCOperator mcmcOperator = schedule.getOperator(op);

            double oldScore = currentScore;
            if (usingFullEvaluation) {
                diagnosticDensities = new HashMap<String, Double>();
                fillDensities(likelihood, diagnosticDensities);
            }

            // assert Profiler.startProfile("Store");

            // The current model is stored here in case the proposal fails
            if (currentModel != null) {
                currentModel.storeModelState();
            }

            // assert Profiler.stopProfile("Store");

            boolean operatorSucceeded = true;
            double hastingsRatio = 1.0;
            boolean accept = false;

            logr[0] = -Double.MAX_VALUE;

            // The new model is proposed
            // assert Profiler.startProfile("Operate");

            if (DEBUG) {
                System.out.println("\n>> Iteration: " + currentState);
                System.out.println("\n&& Operator: " + mcmcOperator.getOperatorName());
            }

            if (mcmcOperator instanceof GeneralOperator) {
                hastingsRatio = ((GeneralOperator) mcmcOperator).operate(likelihood);
            } else {
                hastingsRatio = mcmcOperator.operate();
            }

            // assert Profiler.stopProfile("Operate");
            if (hastingsRatio == Double.NEGATIVE_INFINITY) {
                // Should the evaluation be short-cutted?
                // Previously this was set to false if OperatorFailedException was thrown.
                // Now a -Inf HR is returned.
                operatorSucceeded = false;
            }

            double score = Double.NaN;
            double deviation = Double.NaN;

            //    System.err.print("" + currentState + ": ");
            if (operatorSucceeded) {

                // The new model is proposed
                // assert Profiler.startProfile("Evaluate");

                if (DEBUG) {
                    System.out.println("** Evaluate");
                }

                long elapsedTime = 0;
                long calculationCount = 0;
                if (PROFILE) {
                    elapsedTime = System.currentTimeMillis();
                    if (likelihood instanceof Profileable) {
                        calculationCount = ((Profileable) likelihood).getTotalCalculationCount();
                    }
                }

                // The new model is evaluated
                score = evaluate(likelihood);

                boolean EXTRA_CHECK = true;

                if (EXTRA_CHECK) {
                    String before = likelihood instanceof CompoundLikelihood ?
                            ((CompoundLikelihood) likelihood).getDiagnosis() : "";

                    likelihood.makeDirty();
                    double newScore = evaluate(likelihood);

                    String after = likelihood instanceof CompoundLikelihood ?
                            ((CompoundLikelihood) likelihood).getDiagnosis() : "";

                    if (Math.abs(score - newScore) > evaluationTestThreshold) {

                        System.err.println("Latent error:");
                        System.err.println(before);
                        System.err.println();
                        System.err.println(after);

                        System.exit(-1);
                    }
                }

                if (PROFILE) {
                    long duration = System.currentTimeMillis() - elapsedTime;
                    mcmcOperator.addEvaluationTime(duration);
                    long newCalculationCount = (likelihood instanceof Profileable) ?
                            ((Profileable) likelihood).getTotalCalculationCount() : 1;
                    mcmcOperator.addCalculationCount(newCalculationCount - calculationCount);

                    if (DEBUG) {
                        System.out.println("Time: " + duration);
                    }
                }

                Map<String, Double> diagnosticOperatorDensities = null;
                if (usingFullEvaluation) {
                    diagnosticOperatorDensities = new HashMap<String, Double>();
                    fillDensities(likelihood, diagnosticOperatorDensities);
                }

                if (score == Double.NEGATIVE_INFINITY && mcmcOperator instanceof GibbsOperator) {
                    if (!(mcmcOperator instanceof GibbsIndependentNormalDistributionOperator) && !(mcmcOperator instanceof GibbsIndependentGammaOperator) && !(mcmcOperator instanceof GibbsIndependentCoalescentOperator) && !(mcmcOperator instanceof GibbsIndependentJointNormalGammaOperator)) {
                        if (DEBUG) {
                            if (likelihood instanceof CompoundLikelihood) {
                                CompoundLikelihood cpl = (CompoundLikelihood) likelihood;
                                System.out.println("\nDEBUG CompoundLikelihood:");
                                for (Likelihood l : cpl.getLikelihoods()) {
                                    if (l instanceof CompoundLikelihood) {
                                        CompoundLikelihood l2 = (CompoundLikelihood)l;
                                        for (Likelihood lTwo : l2.getLikelihoods()) {
                                            System.out.println("\t" + lTwo.prettyName() + " = " + lTwo.getLogLikelihood());
                                        }
                                    }
                                    System.out.println(l.prettyName() + " = " + l.getLogLikelihood());
                                }
                                System.out.println("Total likelihood = " + cpl.getLogLikelihood() + "\n");
                            }
                        }
                        Logger.getLogger("error").severe("State " + currentState + ": A Gibbs operator, " + mcmcOperator.getOperatorName() + ", returned a state with zero likelihood.");
                    }
                }

                if (score == Double.POSITIVE_INFINITY ||
                        Double.isNaN(score) ) {
                    if (likelihood instanceof CompoundLikelihood) {
                        Logger.getLogger("error").severe("State "+currentState+": A likelihood returned with a numerical error:\n" +
                                ((CompoundLikelihood)likelihood).getDiagnosis());
                    } else {
                        Logger.getLogger("error").severe("State "+currentState+": A likelihood returned with a numerical error.");
                    }

                    if (DEBUG) {
                        if (likelihood instanceof CompoundLikelihood) {
                            CompoundLikelihood cpl = (CompoundLikelihood) likelihood;
                            System.out.println("\nDEBUG CompoundLikelihood:");
                            for (Likelihood l : cpl.getLikelihoods()) {
                                if (l instanceof CompoundLikelihood) {
                                    CompoundLikelihood l2 = (CompoundLikelihood) l;
                                    for (Likelihood lTwo : l2.getLikelihoods()) {
                                        System.out.println("\t" + lTwo.prettyName() + " = " + lTwo.getLogLikelihood());
                                    }
                                }
                                System.out.println(l.prettyName() + " = " + l.getLogLikelihood());
                            }
                            System.out.println("Total likelihood = " + cpl.getLogLikelihood() + "\n");
                        }
                    }

                    // If the user has chosen to ignore this error then we transform it
                    // to a negative infinity so the state is rejected.
                    score = Double.NEGATIVE_INFINITY;
                }

                if (usingFullEvaluation) {

                    // This is a test that the state was correctly evaluated. The
                    // likelihood of all components of the model are flagged as
                    // needing recalculation, then the full likelihood is calculated
                    // again and compared to the first result. This checks that the
                    // BEAST is aware of all changes that the operator induced.

                    likelihood.makeDirty();
                    final double testScore = evaluate(likelihood);

                    Map<String, Double> densitiesAfter = new HashMap<String, Double>();
                    fillDensities(likelihood, densitiesAfter);

                    if (Math.abs(testScore - score) > evaluationTestThreshold) {
                        StringBuilder sb = new StringBuilder();
                        sb.append("State "+currentState+": State was not correctly calculated after an operator move.\n"
                                + "Likelihood evaluation: " + score
                                + "\nFull Likelihood evaluation: " + testScore
                                + "\n" + "Operator: " + mcmcOperator
                                + " " + mcmcOperator.getOperatorName()
                                + "\n\n");

                        sb.append("Discrepancies:\n");
                        for (String key : diagnosticOperatorDensities.keySet()) {
                            if (Math.abs(diagnosticOperatorDensities.get(key) - densitiesAfter.get(key)) > evaluationTestThreshold) {
                                sb.append(key + ": " + diagnosticOperatorDensities.get(key) +
                                        " -> " + densitiesAfter.get(key) + "\n");
                            }
                        }
                        sb.append("\n");

                        Logger.getLogger("error").severe(sb.toString());
                        fullEvaluationError = true;
                    }
                }

                if (score > bestScore) {
                    bestScore = score;
                    fireBestModel(currentState, currentModel);
                }

                accept = mcmcOperator instanceof GibbsOperator || acceptor.accept(oldScore, score, hastingsRatio, logr);

                deviation = score - oldScore;
            }

            // The new model is accepted or rejected
            if (accept) {
                if (DEBUG) {
                    System.out.println("** Move accepted: new score = " + score
                            + ", old score = " + oldScore);
                }

                mcmcOperator.accept(deviation);
                currentModel.acceptModelState();
                currentScore = score;

            } else {
                if (DEBUG) {
                    System.out.println("** Move rejected: new score = " + score
                            + ", old score = " + oldScore + " (logr = " + logr[0] + ")");
                }

                mcmcOperator.reject();

                // assert Profiler.startProfile("Restore");

                currentModel.restoreModelState();

                if (usingFullEvaluation) {
                    // This is a test that the state is correctly restored. The
                    // restored state is fully evaluated and the likelihood compared with
                    // that before the operation was made.

                    likelihood.makeDirty();
                    final double testScore = evaluate(likelihood);

                    Map<String, Double> densitiesAfter = new HashMap<String, Double>();
                    fillDensities(likelihood, densitiesAfter);

                    if (Math.abs(testScore - oldScore) > evaluationTestThreshold) {
                        StringBuilder sb = new StringBuilder();
                        sb.append("State "+currentState+": State was not correctly restored after reject step.\n"
                                + "Likelihood before: " + oldScore
                                + " Likelihood after: " + testScore
                                + "\n" + "Operator: " + mcmcOperator
                                + " " + mcmcOperator.getOperatorName()
                                + "\n\n");

                        sb.append("Discrepancies:\n");
                        for (String key : diagnosticDensities.keySet()) {
                            if (Math.abs(diagnosticDensities.get(key) - densitiesAfter.get(key)) > evaluationTestThreshold) {
                                sb.append(key + ": " + diagnosticDensities.get(key) +
                                        " -> " + densitiesAfter.get(key) + "\n");
                            }
                        }
                        sb.append("\n");

                        Logger.getLogger("error").severe(sb.toString());
                        fullEvaluationError = true;
                    }
                }
            }
            // assert Profiler.stopProfile("Restore");


            if (isAdapting(mcmcOperator)) {
                adaptAcceptanceProbability((AdaptableMCMCOperator) mcmcOperator, logr[0]);
            }

            if (usingFullEvaluation) {
                if (schedule.getMinimumAcceptAndRejectCount() >= minOperatorCountForFullEvaluation &&
                        currentState >= fullEvaluationCount) {
                    // full evaluation is only switched off when each operator has done a
                    // minimum number of operations (currently 1) and fullEvaluationCount
                    // operations in total.

                    usingFullEvaluation = false;
                    if (fullEvaluationError) {
                        // If there has been an error then stop with an error
                        throw new RuntimeException(
                                "One or more evaluation errors occurred during the test phase of this\n" +
                                        "run. These errors imply critical errors which may produce incorrect\n" +
                                        "results.");
                    }
                }
            }

            fireEndCurrentIteration(currentState);

            currentState += 1;
        }

        currentLength = currentState;

        return currentLength;
    }

    public void terminateChain() {
        fireFinished(currentLength);

        // Profiler.report();
    }

    public Likelihood getLikelihood() {
        return likelihood;
    }

    public Model getModel() {
        return likelihood.getModel();
    }

    public OperatorSchedule getSchedule() {
        return schedule;
    }

    public Acceptor getAcceptor() {
        return acceptor;
    }

    public double getInitialScore() {
        return initialScore;
    }

    public double getBestScore() {
        return bestScore;
    }

    public long getCurrentLength() {
        return currentLength;
    }

    public void setCurrentLength(long currentLength) {
        this.currentLength = currentLength;
    }

    public double getCurrentScore() {
        return currentScore;
    }

    public void pleaseStop() {
        pleaseStop = true;
    }

    public boolean isStopped() {
        return isStopped;
    }

    public double evaluate() {
        return evaluate(likelihood);
    }

    protected double evaluate(Likelihood likelihood) {

        double logPosterior = 0.0;

        final double logLikelihood = likelihood.getLogLikelihood();

        if (Double.isNaN(logLikelihood)) {
            return Double.NEGATIVE_INFINITY;
        }
        // System.err.println("** " + logPosterior + " + " + logLikelihood +
        // " = " + (logPosterior + logLikelihood));
        logPosterior += logLikelihood;

        return logPosterior;
    }

    public boolean isAdapting(MCMCOperator operator) {
        return (isCurrentlyAdapting && operator instanceof AdaptableMCMCOperator);
    }

    /**
     * Updates the proposal parameter, based on the target acceptance
     * probability This method relies on the proposal parameter being a
     * decreasing function of acceptance probability.
     *
     * @param op   The operator
     * @param logr The log likelihood ratio
     */
    public void adaptAcceptanceProbability(AdaptableMCMCOperator op, double logr) {

        final boolean isAdaptable = op.getMode() == AdaptationMode.ADAPTATION_ON || (op.getMode() != AdaptationMode.ADAPTATION_OFF);

        if (DEBUG) {
            System.out.println("adaptAcceptanceProbability " + isAdaptable);
        }

        if (isAdaptable) {
            final double p = op.getAdaptableParameter();

            final double i = schedule.getOptimizationTransform().transform(op.getAdaptationCount() + 2);

            double acceptance;
            if (useSmoothedAcceptanceProbability) {
                acceptance = op.getSmoothedAcceptanceProbability();
            } else {
                acceptance = Math.exp(logr);
            }

            final double target = op.getTargetAcceptanceProbability();

            final double newp = p + ((1.0 / i) * (acceptance - target));

            if (newp > -Double.MAX_VALUE && newp < Double.MAX_VALUE) {
                op.setAdaptableParameter(newp);
                if (DEBUG) {
                    System.out.println("Setting coercable parameter: " + newp + " target: " + target + " current: " + acceptance);
                }
            }
        }
    }

    public void addMarkovChainListener(MarkovChainListener listener) {
        if (listener != null) {
            listeners.add(listener);
        }
    }

    public void removeMarkovChainListener(MarkovChainListener listener) {
        listeners.remove(listener);
    }


    private void fireBestModel(long state, Model bestModel) {

        for (MarkovChainListener listener : listeners) {
            listener.bestState(state, this, bestModel);
        }
    }

    private void fireCurrentModel(long state, Model currentModel) {
        for (MarkovChainListener listener : listeners) {
            listener.currentState(state, this, currentModel);
        }
    }

    private void fireFinished(long chainLength) {

        for (MarkovChainListener listener : listeners) {
            listener.finished(chainLength, this);
        }
    }

    private void fireEndCurrentIteration(long state) {
    }

    private void fillDensities(Likelihood like, Map<String, Double> densities) {
        if (like instanceof CompoundLikelihood) {
            for (Likelihood subLike : ((CompoundLikelihood) like).getLikelihoods()) {
                fillDensities(subLike, densities);
            }
        } else {
            final double logLikelihood = like.getLogLikelihood();
            densities.put(like.prettyName(), logLikelihood);
        }
    }

    private final ArrayList<MarkovChainListener> listeners = new ArrayList<MarkovChainListener>();
}
