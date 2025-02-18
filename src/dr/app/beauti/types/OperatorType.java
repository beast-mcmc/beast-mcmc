/*
 * OperatorType.java
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

package dr.app.beauti.types;

import dr.evomodel.operators.BitFlipInSubstitutionModelOperator;
import dr.inference.operators.RateBitExchangeOperator;
import dr.inferencexml.operators.ScaleOperatorParser;

/**
 * @author Alexei Drummond
 */
public enum OperatorType {

    SCALE("scale"),
    RANDOM_WALK("randomWalk"),
    RANDOM_WALK_ABSORBING("randomWalkAbsorbing"),
    RANDOM_WALK_REFLECTING("randomWalkReflecting"),
    RANDOM_WALK_LOG("randomWalkLog"),
    RANDOM_WALK_LOGIT("randomWalkLogit"),
    INTEGER_RANDOM_WALK("integerRandomWalk"),
    UP_DOWN("upDown"),
    SCALE_ALL(ScaleOperatorParser.SCALE_ALL),
    SCALE_INDEPENDENTLY("scaleIndependently"),
    DELTA_EXCHANGE("deltaExchange"),
    WEIGHTED_DELTA_EXCHANGE("weightedDeltaExchange"),
    INTEGER_DELTA_EXCHANGE("integerDeltaExchange"),
    SWAP("swap"),
    BITFLIP("bitFlip"),
    BITFLIP_IN_SUBST(BitFlipInSubstitutionModelOperator.BIT_FLIP_OPERATOR),// bitFlipInSubstitutionModelOperator
    RATE_BIT_EXCHANGE(RateBitExchangeOperator.OPERATOR_NAME), // rateBitExchangeOperator
    TREE_BIT_MOVE("treeBitMove"),
    SAMPLE_NONACTIVE("sampleNoneActiveOperator"),
    SCALE_WITH_INDICATORS("scaleWithIndicators"),
    UNIFORM("uniform"),
    INTEGER_UNIFORM("integerUniform"),
    SUBTREE_SLIDE("subtreeSlide"),
    FIXED_HEIGHT_SUBTREE_PRUNE_REGRAFT("fixedHeightSPR"),
    SUBTREE_LEAP("subtreeLeap"),
    UNIFORM_SUBTREE_PRUNE_REGRAFT("uniformSPG"),
    ADAPTIVE_MULTIVARIATE("adaptiveMultivariate"),
    NARROW_EXCHANGE("narrowExchange"),
    WIDE_EXCHANGE("wideExchange"),
    EMPIRICAL_TREE_SWAP("empiricalSwap"),
    GMRF_BLOCKUPDATE_OPERATOR("gmrfBlockUpdateOperator"),
    SKY_GRID_BLOCKUPDATE_OPERATOR("gmrfBlockUpdateOperator"),
    SKY_GRID_HMC_OPERATOR("gmrfHMCOperator"),
//    PRECISION_GMRF_OPERATOR("precisionGMRFOperator"),
    WILSON_BALDING("wilsonBalding"),
    RELAXED_CLOCK_HMC_RATE_OPERATOR("relaxedClockHMCRateOperator"),
    RELAXED_CLOCK_HMC_SCALE_OPERATOR("relaxedClockHMCScaleOperator"),
    SHRINKAGE_CLOCK_HMC_OPERATOR("shrinkageClockHMCOperator"),
    SHRINKAGE_CLOCK_GIBBS_OPERATOR("relaxedClockGibbsOperator"),

    //BFT operators
    NODE_HEIGHT_OPERATOR_UNIFORM("nodeHeightOperator-uniform"),
    NODE_HEIGHT_OPERATOR_ROOT("nodeHeightOperator-scaleRoot");

    OperatorType(String displayName) {
        this.displayName = displayName;
    }

    public String toString() {
        return displayName;
    }

    private final String displayName;
}
