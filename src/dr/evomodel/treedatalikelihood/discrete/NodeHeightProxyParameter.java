/*
 * NodeHeightProxyParameter.java
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class NodeHeightProxyParameter extends Parameter.Proxy implements Bounds<Double> {

    private TreeModel tree;
    private TreeParameterModel indexHelper;

    public NodeHeightProxyParameter(String name,
                                    TreeModel tree,
                                    boolean includeRoot) {
        super(name, includeRoot ? tree.getInternalNodeCount() : tree.getInternalNodeCount() - 1);
        this.tree = tree;
        this.indexHelper = new TreeParameterModel(tree,
                new Parameter.Default(includeRoot ? tree.getInternalNodeCount() : tree.getInternalNodeCount() - 1, 0.0),
                includeRoot);
                
    }

    public TreeModel getTree() {
        return tree;
    }

    private int getNodeNumber(int index) {
        return indexHelper.getNodeNumberFromParameterIndex(index + tree.getExternalNodeCount());
    }

    @Override
    public double getParameterValue(int dim) {
        return tree.getNodeHeight(tree.getNode(getNodeNumber(dim)));
    }

    @Override
    public void setParameterValue(int dim, double value) {
        tree.setNodeHeight(tree.getNode(getNodeNumber(dim)), value);
        tree.pushTreeChangedEvent(tree.getNode(getNodeNumber(dim)));
    }

    @Override
    public void setParameterValueQuietly(int dim, double value) {
        tree.setNodeHeightQuietly(tree.getNode(getNodeNumber(dim)), value);
    }
    @Override
    public Bounds<Double> getBounds() {
        return this;
    }
    //bounds implementation
    @Override
    public Double getUpperLimit(int dimension) {
        NodeRef node = tree.getNode(getNodeNumber(dimension));
        if(tree.isRoot(node)){
            return Double.POSITIVE_INFINITY;
        }else{
            return tree.getNodeHeight(tree.getParent(node));
        }
        
    }

    @Override
    public Double getLowerLimit(int dimension) {
        // TODO Auto-generated method stub
        NodeRef node = tree.getNode(getNodeNumber(dimension));
        if (tree.isExternal(node)) {
            return 0.0;
        } else {
            double max = 0.0;
            for(int i = 0; i < tree.getChildCount(node); i++) {
                max = Math.max(max, tree.getNodeHeight(tree.getChild(node,i)));
            }
            return max;
        }
    }

    @Override
    public int getBoundsDimension() {
                return this.getDimension();       
    }
    public String toString() {
        StringBuilder buffer = new StringBuilder(String.valueOf(getParameterValue(0)));
        Bounds bounds = null;

        for (int i = 1; i < getDimension(); i++) {
            buffer.append("\t").append(String.valueOf(getParameterValue(i)));
        }
        return buffer.toString();
    }

    @Override
    public void fireParameterChangedEvent() {
        tree.pushTreeChangedEvent(TreeChangedEvent.create(true, true));
    }

    @Override
    public void setParameterValueNotifyChangedAll(int dim, double value) {
        setParameterValue(dim, value);
    }

    private static final String NODE_HEIGHT_PARAMETER = "nodeHeightProxyParameter";
    public static final String INCLUDE_ROOT = "rootNode";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            boolean includeRoot = xo.getBooleanAttribute(INCLUDE_ROOT);
            TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
            return new NodeHeightProxyParameter(NODE_HEIGHT_PARAMETER, tree, includeRoot);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    AttributeRule.newBooleanRule(INCLUDE_ROOT),
                    new ElementRule(TreeModel.class),
            };
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return NodeHeightProxyParameter.class;
        }

        @Override
        public String getParserName() {
            return NODE_HEIGHT_PARAMETER;
        }
    };
}


