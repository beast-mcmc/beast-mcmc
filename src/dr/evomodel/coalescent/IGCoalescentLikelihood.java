/*
 * CoalescentLikelihood.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.evomodel.coalescent;

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.IntervalList;
import dr.evolution.coalescent.IntervalType;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodelxml.coalescent.CoalescentLikelihoodParser;
import dr.math.Binomial;
import dr.math.LogTricks;

import java.util.List;


/**
 * A likelihood function for the coalescent. Takes a tree and a demographic model.
 *
 * Parts of this class were derived from C++ code provided by Oliver Pybus.
 *
 * @version $Id: NewCoalescentLikelihood.java,v 1.6 2005/05/24 20:25:57 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Luiz Max Carvalho
 */
public final class IGCoalescentLikelihood extends AbstractCoalescentLikelihood implements Units {



	// PUBLIC STUFF
	public IGCoalescentLikelihood(Tree tree,
			TaxonList includeSubtree,
			List<TaxonList> excludeSubtrees,
			double alpha, double beta) throws TreeUtils.MissingTaxonException {

		super(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, tree, includeSubtree, excludeSubtrees);

		this.alpha = alpha;
		this.beta = beta;
	}


	/**
	 * Calculates the log likelihood of this set of coalescent intervals,
	 * given doubles alpha and beta
	 */
	public double calculateLogLikelihood() {
		double lnL =  calculateLogLikelihood(getIntervals(), alpha, beta);
		return lnL;
	}

	public static double calculateLogLikelihood(IntervalList intervals, double alpha, double beta) {

		double logL = 0.0;
		double ldens = 0.0;		
		final int n = intervals.getIntervalCount();
		for (int i = 0; i < n; i++) {

			final double duration = intervals.getInterval(i);
			final int lineageCount = intervals.getLineageCount(i);
			final double kChoose2 = lineageCount*(lineageCount-1.0)/2.0;

			if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {
				ldens = marginalLogDensity(duration, kChoose2, alpha, beta);
				logL +=  ldens;
			}else {
				// sampling interval
				if(lineageCount == 1.0 || duration == 0.0) {
					ldens = 0.0;
					logL += ldens;
				}else {
					ldens = invGammaNoCoalescentlogProb(duration, kChoose2, alpha, beta);
					logL += ldens;
				}
			}
			//		System.err.println("i=" + i + " k=" + lineageCount + " delta= " + duration + " isCoal=" + intervals.getIntervalType(i) +  " ldens= " + ldens + "\n");
		}
//				System.err.println("logLik=" + logL + "\n");
		return logL;
	}

	public static double marginalLogDensity(double t, double coeff, double alpha, double beta) {
		double lnum = Math.log(coeff) + Math.log(alpha) + alpha * Math.log(beta);
		double ldenom = (alpha + 1) * Math.log(coeff * t + beta);
		return(lnum - ldenom);
	}

	public static double invGammaNoCoalescentlogProb(double x, double coeff, double alpha, double beta) {
		//Gives (log) Pr(no coalescence | time = x, coeff, alpha, beta)
		// Pr(no coal | time =x, coeff, Ne) = exp(-x * coeff/Ne). Then integrate against prior on Ne.
		double logL = alpha * ( Math.log(beta)-Math.log(coeff*x + beta) );
		return logL;
	}

	@Override
	public Type getUnits() {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public void setUnits(Type units) {
		// TODO Auto-generated method stub			
	}

	// PRIVATE STUFF 
	private double alpha;
	private double beta;
}