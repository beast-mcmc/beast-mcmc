/*
 * SiteModelParser.java
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

package dr.evomodelxml.siteratemodel;

import java.util.logging.Logger;

import dr.evomodel.siteratemodel.DiscretizedSiteRateModel;
import dr.evomodel.siteratemodel.GammaSiteRateDelegate;
import dr.evomodel.siteratemodel.HomogeneousRateDelegate;
import dr.evomodel.siteratemodel.SiteRateDelegate;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;
import dr.xml.XORRule;

/**
 * This is a replacement to GammaSiteModelParser to keep old XML that used
 * the <siteModel></siteModel> element working.
 * @author Andrew Rambaut
 */
public class SiteModelParser extends AbstractXMLObjectParser {

    public static final String SITE_MODEL = "siteModel";
    public static final String SUBSTITUTION_MODEL = "substitutionModel";
    public static final String MUTATION_RATE = "mutationRate";
    public static final String SUBSTITUTION_RATE = "substitutionRate";
    public static final String RELATIVE_RATE = "relativeRate";
    public static final String WEIGHT = "weight";
    public static final String GAMMA_SHAPE = "gammaShape";
    public static final String GAMMA_CATEGORIES = "gammaCategories";
    public static final String PROPORTION_INVARIANT = "proportionInvariant";

    public String getParserName() {
        return  SITE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String msg = "";
        SubstitutionModel substitutionModel;

        double muWeight = 1.0;

        Parameter muParam = null;
        if (xo.hasChildNamed(SUBSTITUTION_RATE)) {
            muParam = (Parameter) xo.getElementFirstChild(SUBSTITUTION_RATE);

            msg += "\n  with initial substitution rate = " + muParam.getParameterValue(0);
        } else  if (xo.hasChildNamed(MUTATION_RATE)) {
            muParam = (Parameter) xo.getElementFirstChild(MUTATION_RATE);

            msg += "\n  with initial substitution rate = " + muParam.getParameterValue(0);
        } else if (xo.hasChildNamed(RELATIVE_RATE)) {
            XMLObject cxo = xo.getChild(RELATIVE_RATE);
            muParam = (Parameter) cxo.getChild(Parameter.class);
            msg += "\n  with initial relative rate = " + muParam.getParameterValue(0);
            if (cxo.hasAttribute(WEIGHT)) {
                muWeight = cxo.getDoubleAttribute(WEIGHT);
                msg += " with weight: " + muWeight;
            }
        }

        Parameter shapeParam = null;
        int catCount = 4;
        if (xo.hasChildNamed(GAMMA_SHAPE)) {
            XMLObject cxo = xo.getChild(GAMMA_SHAPE);
            catCount = cxo.getIntegerAttribute(GAMMA_CATEGORIES);

            shapeParam = (Parameter) cxo.getChild(Parameter.class);

            msg += "\n  " + catCount + " category discrete gamma with initial shape = " + shapeParam.getParameterValue(0);
            msg += "\n  using equal weight discretization of gamma distribution";
        }

        Parameter invarParam = null;
        if (xo.hasChildNamed(PROPORTION_INVARIANT)) {
            invarParam = (Parameter) xo.getElementFirstChild(PROPORTION_INVARIANT);
            msg += "\n  initial proportion of invariant sites = " + invarParam.getParameterValue(0);
        }

        if (msg.length() > 0) {
            Logger.getLogger("dr.evomodel").info("\nCreating site rate model: " + msg);
        } else {
            Logger.getLogger("dr.evomodel").info("\nCreating site rate model.");
        }

        SiteRateDelegate delegate;
        if (shapeParam != null || invarParam != null) {
            delegate = new GammaSiteRateDelegate("GammaSiteRateDelegate", shapeParam, catCount, invarParam);
        } else {
            delegate = new HomogeneousRateDelegate("HomogeneousRateDelegate");
        }

        DiscretizedSiteRateModel siteRateModel = new DiscretizedSiteRateModel(SiteModel.SITE_MODEL, muParam, muWeight, delegate);

        if (xo.hasChildNamed(SUBSTITUTION_MODEL)) {

//        	System.err.println("Doing the substitution model stuff");

            // set this to pass it along to the OldTreeLikelihoodParser...
            substitutionModel = (SubstitutionModel) xo.getElementFirstChild(SUBSTITUTION_MODEL);
            siteRateModel.setSubstitutionModel(substitutionModel);

        }

        return siteRateModel;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A DiscretizedSiteRateModel that has a gamma distributed rates across sites";
    }

    @Override
    public String[] getParserNames() {
        return super.getParserNames();
    }

    public Class<DiscretizedSiteRateModel> getReturnType() {
        return DiscretizedSiteRateModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {

            new ElementRule(SUBSTITUTION_MODEL, new XMLSyntaxRule[]{
                    new ElementRule(SubstitutionModel.class)
            }, true),

            new XORRule(
                    new XORRule(
                            new ElementRule(SUBSTITUTION_RATE, new XMLSyntaxRule[]{
                                    new ElementRule(Parameter.class)
                            }),
                            new ElementRule(MUTATION_RATE, new XMLSyntaxRule[]{
                                    new ElementRule(Parameter.class)
                            })
                    ),
                    new ElementRule(RELATIVE_RATE, new XMLSyntaxRule[]{
                            AttributeRule.newDoubleRule(WEIGHT, true),
                            new ElementRule(Parameter.class)
                    }), true
            ),

            new ElementRule(GAMMA_SHAPE, new XMLSyntaxRule[]{
                    AttributeRule.newIntegerRule(GAMMA_CATEGORIES, true),
                    new ElementRule(Parameter.class)
            }, true),

            new ElementRule(PROPORTION_INVARIANT, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true)
    };

}//END: class
