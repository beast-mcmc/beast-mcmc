/*
 * NodeHeightOperatorParser.java
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

import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.operators.RandomWalkNodeHeightOperator;
import dr.evomodel.operators.ScaleNodeHeightOperator;
import dr.evomodel.operators.UniformNodeHeightOperator;
import dr.evomodel.tree.TreeModel;
//import dr.evomodel.treelikelihood.thorneytreelikelihood.MultiMoveUniformNodeHeightOperator;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 * @author Andrew Rambaut
 */
public class NodeHeightOperatorParser extends AbstractXMLObjectParser {

    public enum OperatorType {
        UNIFORM("uniform"),
        RANDOMWALK("random walk"),
        SCALEROOT("scale root"),
        SCALEALL("scale all internal");


        OperatorType(String name) {
            this.name = name;
        }

        @Override
        public String toString() {
            return name;
        }

        String name;
    }

    public static final String NODE_HEIGHT_OPERATOR = "nodeHeightOperator";

    public static final String SIZE = "size";
    public static final String SCALE_FACTOR = "scaleFactor";
    public static final String MEAN_COUNT = "meanCount";
    public static final String TARGET_ACCEPTANCE = "targetAcceptance";
    public static final String OPERATOR_TYPE = "type";

    public String getParserName() {
        return NODE_HEIGHT_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        OperatorType operatorType = OperatorType.UNIFORM;
        if (xo.hasAttribute(OPERATOR_TYPE)) {
            try {
                operatorType = OperatorType.valueOf(xo.getStringAttribute(OPERATOR_TYPE).trim().toUpperCase());
            } catch (IllegalArgumentException iae) {
                throw new XMLParseException("Unrecognised operator type attribute: " + xo.getStringAttribute(OPERATOR_TYPE));
            }
        }

        if (operatorType != OperatorType.UNIFORM) {
            AdaptationMode mode = AdaptationMode.parseMode(xo);

            double tuningParameter = 0.75;

            if (xo.hasAttribute(SIZE)) {
                tuningParameter = xo.getDoubleAttribute(SIZE);
                if (tuningParameter <= 0.0) {
                    throw new XMLParseException("The UniformNodeHeightOperator size attribute must be positive and non-zero.");
                }
            }
            if (xo.hasAttribute(SCALE_FACTOR)) {
                tuningParameter = xo.getDoubleAttribute(SCALE_FACTOR);
                if (tuningParameter <= 0.0 || tuningParameter >= 1.0) {
                    throw new XMLParseException("The UniformNodeHeightOperator scaleFactor attribute must be between 0 and 1.");
                }
            }

            final double targetAcceptance = xo.getAttribute(TARGET_ACCEPTANCE, 0.234);

            if (targetAcceptance <= 0.0 || targetAcceptance >= 1.0) {
                throw new XMLParseException("Target acceptance probability has to lie in (0, 1)");
            }

            switch (operatorType) {
                case RANDOMWALK:
                    return new RandomWalkNodeHeightOperator(treeModel, weight, tuningParameter, mode, targetAcceptance);
                case SCALEROOT:
                case SCALEALL:
                    return new ScaleNodeHeightOperator(treeModel, weight, tuningParameter, operatorType, mode, targetAcceptance);
                default:
                    throw new IllegalArgumentException("Unknown operator type");
            }
        }

        if (xo.hasAttribute(SIZE)|| xo.hasAttribute(SCALE_FACTOR) || xo.hasAttribute(TARGET_ACCEPTANCE) || xo.hasAttribute(AdaptableMCMCOperator.AUTO_OPTIMIZE)) {
            throw new XMLParseException("Uniform node height operator does not take adaptable operator arguments");
        }

        return new UniformNodeHeightOperator(treeModel, weight);
    }

    public String getParserDescription() {
        return "An operator that moves node heights about.";
    }

    public Class getReturnType() {
        return AbstractTreeOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newDoubleRule(SIZE, true),
            AttributeRule.newDoubleRule(SCALE_FACTOR, true),
            AttributeRule.newDoubleRule(TARGET_ACCEPTANCE, true),
            AttributeRule.newStringRule(OPERATOR_TYPE, true),
            AttributeRule.newBooleanRule(AdaptableMCMCOperator.AUTO_OPTIMIZE, true),
            new ElementRule(TreeModel.class)
    };

}