/*
 * IGCoalescentLikelihoodParser.java
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

package dr.evomodelxml.coalescent;

import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.coalescent.IGCoalescentLikelihood;
import dr.evomodel.coalescent.DemographicModel;
import dr.evomodel.coalescent.MultiLociTreeSet;
import dr.evomodel.coalescent.OldAbstractCoalescentLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 */
public class IGCoalescentLikelihoodParser extends AbstractXMLObjectParser {

    public static final String IGCOALESCENT_LIKELIHOOD = "coalescentLikelihoodIG";
    public static final String MODEL = "model";
    public static final String POPULATION_TREE = "populationTree";
    public static final String POPULATION_FACTOR = "factor";
    
    public static final String ALPHA = "alpha";
    public static final String BETA = "beta";

    public static final String INCLUDE = "include";
    public static final String EXCLUDE = "exclude";

    public String getParserName() {
        return IGCOALESCENT_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(MODEL);
        DemographicModel demoModel = null;
        
        double alpha = xo.getDoubleAttribute(ALPHA);
        double beta = xo.getDoubleAttribute(BETA);

        List<TreeModel> trees = new ArrayList<TreeModel>();
        List<Double> popFactors = new ArrayList<Double>();
        MultiLociTreeSet treesSet = demoModel instanceof MultiLociTreeSet ? (MultiLociTreeSet) demoModel : null;

        for (int k = 0; k < xo.getChildCount(); ++k) {
            final Object child = xo.getChild(k);
            if (child instanceof XMLObject) {
                cxo = (XMLObject) child;
                if (cxo.getName().equals(POPULATION_TREE)) {
                    final TreeModel t = (TreeModel) cxo.getChild(TreeModel.class);
                    assert t != null;
                    trees.add(t);

                    popFactors.add(cxo.getAttribute(POPULATION_FACTOR, 1.0));
                }
            }
//                in the future we may have arbitrary multi-loci element
//                else if( child instanceof MultiLociTreeSet )  {
//                    treesSet = (MultiLociTreeSet)child;
//                }
        }

        TreeModel treeModel = null;
        if (trees.size() == 1 && popFactors.get(0) == 1.0) {
            treeModel = trees.get(0);
        } else if (trees.size() > 1) {
            treesSet = new MultiLociTreeSet.Default(trees, popFactors);
        } else if (!(trees.size() == 0 && treesSet != null)) {
            throw new XMLParseException("Incorrectly constructed likelihood element");
        }

        TaxonList includeSubtree = null;

        if (xo.hasChildNamed(INCLUDE)) {
            includeSubtree = (TaxonList) xo.getElementFirstChild(INCLUDE);
        }

        List<TaxonList> excludeSubtrees = new ArrayList<TaxonList>();

        if (xo.hasChildNamed(EXCLUDE)) {
            cxo = xo.getChild(EXCLUDE);
            for (int i = 0; i < cxo.getChildCount(); i++) {
                excludeSubtrees.add((TaxonList) cxo.getChild(i));
            }
        }

        if (treeModel != null) {
            try {
                return new IGCoalescentLikelihood(treeModel, includeSubtree, excludeSubtrees, alpha, beta);
            } catch (TreeUtils.MissingTaxonException mte) {
                throw new XMLParseException("treeModel missing a taxon from taxon list in " + getParserName() + " element");
            }
        } else {
            if (includeSubtree != null || excludeSubtrees.size() > 0) {
                throw new XMLParseException("Include/Exclude taxa not supported for multi locus sets");
            }
            // Use old code for multi locus sets.
            // This is a little unfortunate but the current code is using AbstractCoalescentLikelihood as
            // a base - and modifing it will probsbly result in a bigger mess.
            return new OldAbstractCoalescentLikelihood(treesSet, demoModel);
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of the tree given the demographic function.";
    }

    public Class getReturnType() {
        return IGCoalescentLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MODEL, new XMLSyntaxRule[]{
                    new ElementRule(DemographicModel.class)
            }, "The demographic model which describes the effective population size over time"),

            new ElementRule(POPULATION_TREE, new XMLSyntaxRule[]{
                    AttributeRule.newDoubleRule(POPULATION_FACTOR, true),
                    new ElementRule(TreeModel.class)
            }, "Tree(s) to compute likelihood for", 0, Integer.MAX_VALUE),

            new ElementRule(INCLUDE, new XMLSyntaxRule[]{
                    new ElementRule(Taxa.class)
            }, "An optional subset of taxa on which to calculate the likelihood (should be monophyletic)", true),

            new ElementRule(EXCLUDE, new XMLSyntaxRule[]{
                    new ElementRule(Taxa.class, 1, Integer.MAX_VALUE)
            }, "One or more subsets of taxa which should be excluded from calculate the likelihood (should be monophyletic)", true)
    };
}
