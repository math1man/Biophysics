package com.ariweiland.biophysics.sampler;

import com.ariweiland.biophysics.lattice.Lattice;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Ari Weiland
 */
public abstract class WangLandauSampler extends Sampler {

    public static final MathContext MC = MathContext.DECIMAL64;
    public static final double F_FINAL = 0.00000001; // 10^-8
    public static final int MIN_H = 20;

    private double flatness = 0.8;  // must be between 0 and 1 exclusive

    protected final Map<Double, BigDecimal> g = new HashMap<>();
    protected final Map<Double, Integer> h = new HashMap<>();

    public WangLandauSampler() {}

    public WangLandauSampler(double flatness) {
        this.flatness = flatness;
    }

    public double getFlatness() {
        return flatness;
    }

    public void setFlatness(double flatness) {
        this.flatness = flatness;
    }

    protected boolean isSufficientlyFlat() {
        if (g.isEmpty() || h.size() != g.size()) {
            return false;
        }
        double total = 0;
        for (int i : h.values()) {
            if (i < MIN_H) {
                return false;
            }
            total += i;
        }
        double average = total / h.size();
        for (int i : h.values()) {
            if (i < flatness * average) {
                return false;
            }
        }
        return true;
    }

    protected void updateMaps(double energy, BigDecimal f) {
        g.put(energy, g(energy).multiply(f, MC));
        if (h.containsKey(energy)) {
            h.put(energy, h.get(energy) + 1);
        } else {
            h.put(energy, 1);
        }
    }

    protected BigDecimal g(double energy) {
        if (!g.containsKey(energy) || g.get(energy).equals(BigDecimal.ZERO)) {
            BigDecimal min = BigDecimal.ZERO;
            for (BigDecimal d : g.values()) {
                if (min.equals(BigDecimal.ZERO) || d.compareTo(BigDecimal.ZERO) > 0 && d.compareTo(min) < 0) {
                    min = d;
                }
            }
            if (min.equals(BigDecimal.ZERO)) {
                min = BigDecimal.ONE; // the starting value
            }
            g.put(energy, min);
        }
        return g.get(energy);
    }

    protected double threshold(Lattice old, Lattice trial, double detailedBalance) {
        return g(old.getEnergy()).divide(g(trial.getEnergy()), MC).doubleValue() * detailedBalance;
    }

    protected void reduceG() {
        BigDecimal min = g.get(0.0);
        for (BigDecimal d : g.values()) {
            if (d.compareTo(min) < 0) {
                min = d;
            }
        }
        for (double e : g.keySet()) {
            g.put(e, g.get(e).divide(min, MC));
        }
    }

    protected Map<Double, Double> convertG() {
        reduceG();
        Map<Double, Double> output = new HashMap<>();
        for (double e : g.keySet()) {
            output.put(e, g.get(e).doubleValue());
        }
        return output;
    }

    protected static class FValue {
        private double f;
        private BigDecimal bdf;

        public FValue(double f) {
            this.f = f;
            this.bdf = BigDecimal.valueOf(f);
        }

        public double asDouble() {
            return f;
        }

        public BigDecimal asBigDecimal() {
            return bdf;
        }

        public void sqrt() {
            f = Math.sqrt(f);
            bdf = BigDecimal.valueOf(f);
        }

        @Override
        public String toString() {
            return String.valueOf(f);
        }
    }
}
