package com.ariweiland.biophysics.research;

import java.util.Map;

/**
 * @author Ari Weiland
 */
public abstract class ErrorFunction {

    public abstract double error(Map<Double, Double> expected, Map<Double, Double> outcome);

    public static final ErrorFunction MEAN_SQUARED = new ErrorFunction() {
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
                if (e < 1) {
                    double o = outcome.get(energy);
                    double ratio = Math.log(o / e);
                    sum += ratio * ratio / -Math.log(e);
                }
            }
            return sum;
        }
    };
}
