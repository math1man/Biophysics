package com.ariweiland.biophysics;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Ari Weiland
 */
public class SpecificHeat {

    public static double boltzmannFactor(double g, double e, double t) {
        return g * Math.exp(-e/t);
    }

    public static double partitionFunction(Map<Double, Double> density, double t) {
        double sum = 0;
        for (double e : density.keySet()) {
            sum += boltzmannFactor(density.get(e), e, t);
        }
        return sum;
    }

    public static double averageEnergy(Map<Double, Double> density, double t) {
        double sum = 0;
        for (double e : density.keySet()) {
            sum += e * boltzmannFactor(density.get(e), e, t);
        }
        return sum / partitionFunction(density, t);
    }

    public static double averageSquaredEnergy(Map<Double, Double> density, double t) {
        double sum = 0;
        for (double e : density.keySet()) {
            sum += e * e * boltzmannFactor(density.get(e), e, t);
        }
        return sum / partitionFunction(density, t);
    }

    public static double specificHeat(Map<Double, Double> density, double t, int n) {
        double avgEnergy = averageEnergy(density, t);
        double avgSquaredEnergy = averageSquaredEnergy(density, t);
        return (avgSquaredEnergy - avgEnergy * avgEnergy) / t / t / n;
    }

    public static Map<Double, Double> specificHeatGraph(Map<Double, Double> density, double maxT, int n) {
        Map<Double, Double> specificHeatGraph = new HashMap<>();
        for (int i = 0; i < 1000; i++) {
            double t = (i + 1) * maxT / 1000;
            double sh = SpecificHeat.specificHeat(density, t, n);
            if (Double.isFinite(sh)) {
                specificHeatGraph.put(t, sh);
            }
        }
        return specificHeatGraph;
    }

}
