/*
 * SubtreeSlideOperatorParser.java
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

package dr.evomodelxml.operators;

import dr.evomodel.bigfasttree.BigFastTreeModel;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class SubtreeSlideOperatorParser extends AbstractXMLObjectParser {

    public static final String SUBTREE_SLIDE = "subtreeSlide";
    public static final String SWAP_RATES = "swapInRandomRate";
    public static final String SWAP_TRAITS = "swapInRandomTrait";
    public static final String DIRICHLET_BRANCHES = "branchesAreScaledDirichlet";
    public static final String TARGET_ACCEPTANCE = "targetAcceptance";

    public static final String TRAIT = "trait";

    public String getParserName() {
        return SUBTREE_SLIDE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean swapRates = xo.getAttribute(SWAP_RATES, false);
        boolean swapTraits = xo.getAttribute(SWAP_TRAITS, false);
        boolean scaledDirichletBranches = xo.getAttribute(DIRICHLET_BRANCHES, false);

        AdaptationMode mode = AdaptationMode.DEFAULT;
        if (xo.hasAttribute(AdaptableMCMCOperator.AUTO_OPTIMIZE)) {
            if (xo.getBooleanAttribute(AdaptableMCMCOperator.AUTO_OPTIMIZE)) {
                mode = AdaptationMode.ADAPTATION_ON;
            } else {
                mode = AdaptationMode.ADAPTATION_OFF;
            }
        }

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

        if((swapRates || swapTraits) && treeModel instanceof BigFastTreeModel) {
            throw new XMLParseException("Swapping rates or traits is not implemented for BFT model yet.");
        }
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        final double targetAcceptance = xo.getAttribute(TARGET_ACCEPTANCE, 0.234);

        final double size = xo.getAttribute("size", 1.0);

        if (Double.isInfinite(size) || size <= 0.0) {
            throw new XMLParseException("size attribute must be positive and not infinite. was " + size +
           " for tree " + treeModel.getId() );
        }

        final boolean gaussian = xo.getBooleanAttribute("gaussian");
        SubtreeSlideOperator operator = new SubtreeSlideOperator(treeModel, weight, size, gaussian,
                swapRates, swapTraits, scaledDirichletBranches, mode, targetAcceptance);

        return operator;
    }

    public String getParserDescription() {
        return "An operator that slides a subtree.";
    }

    public Class getReturnType() {
        return SubtreeSlideOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            // Make size optional. If not given or equals zero, size is set to half of average tree branch length.
            AttributeRule.newDoubleRule("size", true),
            AttributeRule.newDoubleRule(TARGET_ACCEPTANCE, true),
            AttributeRule.newBooleanRule("gaussian"),
            AttributeRule.newBooleanRule(SWAP_RATES, true),
            AttributeRule.newBooleanRule(SWAP_TRAITS, true),
            AttributeRule.newBooleanRule(AdaptableMCMCOperator.AUTO_OPTIMIZE, true),
            new ElementRule(TreeModel.class)
    };

}
