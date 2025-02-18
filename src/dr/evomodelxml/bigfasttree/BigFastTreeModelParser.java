/*
 * BigFastTreeModelParser.java
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

package dr.evomodelxml.bigfasttree;

import dr.evolution.tree.Tree;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evomodel.bigfasttree.BigFastTreeModel;
import dr.evomodel.bigfasttree.thorney.RootHeightProxyParameter;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightProxyParameter;
import dr.evomodelxml.tree.TreeModelParser;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterParser;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 */
public class BigFastTreeModelParser extends AbstractXMLObjectParser {


    public BigFastTreeModelParser() {
        rules = new XMLSyntaxRule[]{
                new ElementRule(Tree.class),
                new ElementRule(TreeModelParser.ROOT_HEIGHT, Parameter.class, "A parameter definition with id only (cannot be a reference!)", false),
                new ElementRule(TreeModelParser.NODE_HEIGHTS,
                            new XMLSyntaxRule[]{
                                    AttributeRule.newBooleanRule(TreeModelParser.ROOT_NODE, true, "If true the root height is included in the parameter"),
                                    AttributeRule.newBooleanRule(TreeModelParser.INTERNAL_NODES, true, "If true the internal node heights (minus the root) are included in the parameter"),
                                    new ElementRule(Parameter.class, "A parameter definition with id only (cannot be a reference!)")
                            }, 1, Integer.MAX_VALUE),
        };
    }

    public String getParserName() {
        return BigFastTreeModel.BIG_FAST_TREE_MODEL;
    }

    /**
     * @return a tree object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Tree tree = (Tree) xo.getChild(Tree.class);
        boolean fixHeights = xo.getAttribute(TreeModelParser.FIX_HEIGHTS, false);
        boolean fixTree = xo.getAttribute(TreeModelParser.FIX_TREE, false);

        BigFastTreeModel treeModel = new BigFastTreeModel(xo.getId(), tree, fixHeights, fixTree);

        // Make proxy parameters
        RootHeightProxyParameter rootHeightProxyParameter = new RootHeightProxyParameter("Placeholder_Root_Proxy",
                treeModel); // id overwritten below
        NodeHeightProxyParameter nodeHeightProxyParameter = new NodeHeightProxyParameter("Placeholder_nodeHeight_Proxy",
                treeModel, false);
        NodeHeightProxyParameter allNodeHeightProxyParameter = new NodeHeightProxyParameter(
                "Placeholder_allNodeHeight_Proxy", treeModel, true);
        // parse proxy parameters
        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof XMLObject) {

                XMLObject cxo = (XMLObject) xo.getChild(i);

                if (cxo.getName().equals(TreeModelParser.ROOT_HEIGHT)) {
                    ParameterParser.replaceParameter(cxo, rootHeightProxyParameter);
                } else if (cxo.getName().equals(TreeModelParser.NODE_HEIGHTS)) {

                    boolean rootNode = cxo.getAttribute(TreeModelParser.ROOT_NODE, false);
                    boolean internalNodes = cxo.getAttribute(TreeModelParser.INTERNAL_NODES, false);
                    boolean leafNodes = cxo.getAttribute(TreeModelParser.LEAF_NODES, false);
                    if (leafNodes) {
                        throw new XMLParseException("leafNodes not currently implemented for bft tree");
                    }
                    if (!rootNode && !internalNodes) {
                        throw new XMLParseException(
                                "internal or internal and rootNode nodes must be selected for the nodeHeights element");
                    }
                    if (internalNodes && rootNode) {
                        ParameterParser.replaceParameter(cxo, allNodeHeightProxyParameter);
                    } else if (internalNodes && !rootNode) {
                        ParameterParser.replaceParameter(cxo, nodeHeightProxyParameter);
                    }
// options for default tree that are not implemented here yet.
                }else if(cxo.getName().equals(TreeModelParser.LEAF_HEIGHT) || cxo.getName().equals(TreeModelParser.LEAF_HEIGHTS)) {
                    throw new XMLParseException("Leaf sampling not implemented for bft trees");
                }else if(cxo.getName().equals(TreeModelParser.NODE_RATES) || cxo.getName().equals(TreeModelParser.NODE_TRAITS)|| cxo.getName().equals(TreeModelParser.LEAF_TRAIT)) {
                    throw new XMLParseException("Leaf sampling not implemented for bft trees");
                }
            }
        }

        Logger.getLogger("dr.evomodel").info("\nCreating the big fast tree model, '" + xo.getId() + "'");


        Logger.getLogger("dr.evomodel").info("  taxon count = " + treeModel.getExternalNodeCount());
        Logger.getLogger("dr.evomodel").info("  tree height = " + treeModel.getNodeHeight(treeModel.getRoot()));
        return treeModel;
    }

    private void setUncertaintyBounds(Parameter newParameter, Taxon taxon) {
        Date date = taxon.getDate();
        if (date != null) {
            double uncertainty = date.getUncertainty();
            if (uncertainty > 0.0) {
                // taxon date not specified to exact value so add appropriate bounds
                double upper = Taxon.getHeightFromDate(date);
                double lower = Taxon.getHeightFromDate(date);
                if (date.isBackwards()) {
                    upper += uncertainty;
                } else {
                    lower -= uncertainty;
                }

                // set the bounds for the given precision
                newParameter.addBounds(new Parameter.DefaultBounds(upper, lower, 1));

                // set the initial value to be mid-point
                newParameter.setParameterValue(0, (upper + lower) / 2);
            }
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents a model of the tree. The tree model includes and attributes of the nodes " +
                "including the age (or <i>height</i>) and the rate of evolution at each node in the tree.";
    }

    public String getExample() {
        return
//                "<!-- the tree model as special sockets for attaching parameters to various aspects of the tree     -->\n" +
//                        "<!-- The treeModel below shows the standard setup with a parameter associated with the root height -->\n" +
//                        "<!-- a parameter associated with the internal node heights (minus the root height) and             -->\n" +
//                        "<!-- a parameter associates with all the internal node heights                                     -->\n" +
//                        "<!-- Notice that these parameters are overlapping                                                  -->\n" +
//                        "<!-- The parameters are subsequently used in operators to propose changes to the tree node heights -->\n" +
                        "<treeModel id=\"treeModel1\">\n" +
                        "	<tree idref=\"startingTree\"/>\n" +
//                        "	<rootHeight>\n" +
//                        "		<parameter id=\"treeModel1.rootHeight\"/>\n" +
//                        "	</rootHeight>\n" +
//                        "	<nodeHeights internalNodes=\"true\" rootNode=\"false\">\n" +
//                        "		<parameter id=\"treeModel1.internalNodeHeights\"/>\n" +
//                        "	</nodeHeights>\n" +
//                        "	<nodeHeights internalNodes=\"true\" rootNode=\"true\">\n" +
//                        "		<parameter id=\"treeModel1.allInternalNodeHeights\"/>\n" +
//                        "	</nodeHeights>\n" +
                        "</treeModel>";

    }

    public Class getReturnType() {
        return BigFastTreeModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules;
}
