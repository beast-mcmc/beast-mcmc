/*
 * InverseGammaDistribution.java
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

package dr.math.distributions;

import dr.math.GammaFunction;
import dr.math.UnivariateFunction;

/**
 * inverse gamma distribution.
 * <p/>
 * (Parameters: shape, scale; mean: ??; variance: ??)
 *
 * @author Joseph Heled
 */
public class InverseGammaDistribution implements Distribution {

    private double shape, scale;

    private final double factor;
    private final double logFactor;

    public InverseGammaDistribution(double shape, double scale) {
        this.shape = shape;
        this.scale = scale;
        this.factor = Math.pow(scale, shape) / Math.exp(GammaFunction.lnGamma(shape));
        this.logFactor = shape * Math.log(scale) - GammaFunction.lnGamma(shape);
    }

    public double getShape() {
        return shape;
    }

    public void setShape(double value) {
        shape = value;
    }

    public double getScale() {
        return scale;
    }

    public void setScale(double value) {
        scale = value;
    }

    public double pdf(double x) {
        return pdf(x, shape, scale, factor);
    }

    public double logPdf(double x) {
        return logPdf(x, shape, scale, logFactor);
    }

    public double cdf(double x) {
        return cdf(x, shape, scale);
    }

    public double quantile(double y) {
        return quantile(y, shape, scale);
    }

    public double mean() {
        return mean(shape, scale);
    }

    public double variance() {
        return variance(shape, scale);
    }

    public double nextInverseGamma() {
        return nextInverseGamma(shape, scale);
    }

    public final UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private final UnivariateFunction pdfFunction = new UnivariateFunction() {
        public final double evaluate(double x) {
            return pdf(x);
        }

        public final double getLowerBound() {
            return 0.0;
        }

        public final double getUpperBound() {
            return Double.POSITIVE_INFINITY;
        }
    };

    /**
     * probability density function of the Gamma distribution
     *
     * @param x     argument
     * @param shape shape parameter
     * @param scale scale parameter
     * @return pdf value
     */
    public static double pdf(double x, double shape, double scale, double factor) {
        if (x <= 0)
            return 0.0;

        final double a = Math.exp(shape*Math.log(scale) - (shape + 1)*Math.log(x) - (scale/x)
                - GammaFunction.lnGamma(shape));

        return factor * a;
    }

    /**
     * the natural log of the probability density function of the distribution
     *
     * @param x     argument
     * @param shape shape parameter
     * @param scale scale parameter
     * @return log pdf value
     */
    public static double logPdf(double x, double shape, double scale, double lfactor) {
        if (x <= 0)
            return Double.NEGATIVE_INFINITY;

        return  lfactor + shape*Math.log(scale) - (shape + 1)*Math.log(x) - (scale/x) - GammaFunction.lnGamma(shape);
    }

    /**
     * cumulative density function of the Gamma distribution
     *
     * @param x     argument
     * @param shape shape parameter
     * @param scale scale parameter
     * @return cdf value
     */
    public static double cdf(double x, double shape, double scale) {        
        if (x <= 0.0 || shape <= 0.0) {
            return 0.0;
        }
        return GammaFunction.incompleteGammaQ(shape, scale/x);
    }

    /**
     * quantile (inverse cumulative density function) of the Gamma distribution
     *
     * @param y     argument
     * @param shape shape parameter
     * @param scale scale parameter
     * @return icdf value
     */
    public static double quantile(double y, double shape, double scale) {
        // this is what R thinks
        final GammaDistribution g = new GammaDistribution(shape, scale);
        return 1/g.quantile(1-y);
    }

    /**
     * mean of the Gamma distribution
     *
     * @param shape shape parameter
     * @param scale scale parameter
     * @return mean
     */
    public static double mean(double shape, double scale) {
        if( shape > 1 ) {
            return scale / (shape - 1);
        }
        return Double.POSITIVE_INFINITY;
    }

    /**
     * variance of the Gamma distribution.
     *
     * @param shape shape parameter
     * @param scale scale parameter
     * @return variance
     */
    public static double variance(double shape, double scale) {
         if( shape > 2 ) {
            return scale*scale / ((shape - 1)*(scale-1)*(scale-2));
        }
         return Double.POSITIVE_INFINITY;
    }

    /**
     * sample from the Gamma distribution. This could be calculated using
     * quantile, but current algorithm is faster.
     *
     * @param shape shape parameter
     * @param scale scale parameter
     * @return sample
     */
    public static double nextInverseGamma(double shape, double scale) {
        return 1.0 / new GammaDistribution(shape, 1/scale).nextGamma();
    }
}