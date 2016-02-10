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

    public static final MathContext MC = MathContext.DECIMAL32;
    public static final double F_FINAL = 0.00000001; // 10^-8

    private double minEnergy = 0;   // must be less than or equal to 0
    private double flatness = 0.8;  // must be between 0 and 1 exclusive

    protected final Map<Double, BigDecimal> g = new HashMap<>();
    protected final Map<Double, Integer> h = new HashMap<>();

    public WangLandauSampler() {}

    public WangLandauSampler(double minEnergy) {
        this.minEnergy = minEnergy;
    }

    public WangLandauSampler(double minEnergy, double flatness) {
        this.minEnergy = minEnergy;
        this.flatness = flatness;
    }

    public double getMinEnergy() {
        return minEnergy;
    }

    public void setMinEnergy(double minEnergy) {
        this.minEnergy = minEnergy;
    }

    public double getFlatness() {
        return flatness;
    }

    public void setFlatness(double flatness) {
        this.flatness = flatness;
    }

    protected void initializeG() {
        g.clear();
        if (minEnergy != 0) { // if there is a minimum energy specified, initialize g
            for (double e = minEnergy; e <= 0; e++) {
                g.put(e, BigDecimal.ZERO);
            }
        }
    }

    protected boolean isSufficientlyFlat() {
        if (g.isEmpty() || h.size() != g.size()) {
            return false;
        }
        double total = 0;
        for (int i : h.values()) {
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

    protected double calculateThreshold(Lattice old, Lattice trial, double fineDetail) {
        return g(old.getEnergy()).divide(g(trial.getEnergy()), MC)
                .multiply(BigDecimal.valueOf(fineDetail), MC).doubleValue();
    }

    protected Map<Double, Double> convertG() {
        Map<Double, Double> output = new HashMap<>();
        BigDecimal min = g.get(0.0);
        for (BigDecimal bd : g.values()) {
            if (bd.compareTo(min) < 0) {
                min = bd;
            }
        }
        for (double e : g.keySet()) {
            BigDecimal d = g.get(e);
            output.put(e, d.divide(min, MC).doubleValue());
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
    }
}
