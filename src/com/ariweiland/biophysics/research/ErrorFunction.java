package com.ariweiland.biophysics.research;

import java.util.Map;

/**
 * @author Ari Weiland
 */
public abstract class ErrorFunction {

    public abstract double error(Map<Double, Double> expected, Map<Double, Double> outcome);

    public static final ErrorFunction MEAN_SQUARED_LOG = new ErrorFunction() {
        @Override
        public double error(Map<Double, Double> expected, Map<Double, Double> outcome) {
            double sum = 0;
            for (double energy : expected.keySet()) {
                double e = expected.get(energy);
                double o = outcome.get(energy);
                double ratio = Math.log(o / e);
                sum += ratio * ratio;
            }
            return sum / expected.size();

        }
    };

    public static final ErrorFunction CHI_SQUARED = new ErrorFunction() {
        @Override
        public double error(Map<Double, Double> expected, Map<Double, Double> outcome) {
            double sum = 0;
            for (double energy : expected.keySet()) {
                double e = expected.get(energy);
                double o = outcome.get(energy);
                double residual = o - e;
                sum += residual * residual / e; // still biased towards higher probabilities
            }
            return sum / expected.size();
        }
    };

    public static final ErrorFunction ANDERSON_DARLING = new ErrorFunction() {
        @Override
        public double error(Map<Double, Double> expected, Map<Double, Double> outcome) {
            double sum = 0;
            for (double energy : expected.keySet()) {
                double e = expected.get(energy);
                double o = outcome.get(energy);
                double residual = o - e;
                sum += residual * residual / (e * (1 - e)); // distinct from Chi Squared here, adds a relative bias towards low e
                // note that in general the weighting is equal for e1 and e2 where e1 + e2 = 1
            }
            return sum / expected.size();
        }
    };

    public static final ErrorFunction CHI_SQUARED_LOG = new ErrorFunction() {
        @Override
        public double error(Map<Double, Double> expected, Map<Double, Double> outcome) {
            double sum = 0;
            for (double energy : expected.keySet()) {
                double e = expected.get(energy);
                double o = outcome.get(energy);
                double residual = Math.log(o / e);
                sum += residual * residual / -Math.log(e);
            }
            return sum / expected.size();
        }
    };

    public static final ErrorFunction BHATTACHARYYA = new ErrorFunction() {
        @Override
        public double error(Map<Double, Double> expected, Map<Double, Double> outcome) {
            double bc = 0;
            for (double energy : expected.keySet()) {
                double e = expected.get(energy);
                double o = outcome.get(energy);
                bc += Math.sqrt(e * o);
            }
            return -Math.log(bc);
        }
    };

    public static final ErrorFunction KULLBACK_LEIBLER = new ErrorFunction() {
        @Override
        public double error(Map<Double, Double> expected, Map<Double, Double> outcome) {
            double sum = 0;
            for (double energy : expected.keySet()) {
                double e = expected.get(energy);
                double o = outcome.get(energy);
                double ratio = Math.log(o / e);
                sum += ratio * e; // biased towards errors in higher probabilities
            }
            return sum;
        }
    };

}
