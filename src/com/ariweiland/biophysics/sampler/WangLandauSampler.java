package com.ariweiland.biophysics.sampler;

import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.lattice.Lattice;
import com.ariweiland.biophysics.peptide.Polypeptide;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Ari Weiland
 */
public class WangLandauSampler extends Sampler {

    public static final int MIN_SAMPLES = 100;
    public static final double FLATNESS = 0.8; // must be between 0 and 1 exclusive
    public static final double F_FINAL = Math.exp(0.00000001); // e^(10^-8)

    private final Map<Double, Double> g = new HashMap<>();
    private final Map<Double, Integer> h = new HashMap<>();

    private boolean accept(Lattice old, Lattice trial) {
        double threshold = g.get(old.getEnergy()) / g.get(trial.getEnergy());
        // potentially slightly faster because randoms do not always need to be generated
        return threshold >= 1 || Math.random() < threshold;
    }

    private boolean isSufficientlyFlat() {
        if (h.isEmpty()) {
            return false;
        }
        double total = 0;
        for (int i : h.values()) {
            if (i < MIN_SAMPLES) { // needs to at least have a base number of samples
                return false;
            }
            total += i;
        }
        double average = total / h.size();
        for (int i : h.values()) {
            if (i < FLATNESS * average) {
                return false;
            }
        }
        return true;
    }

    private void updateMaps(double e, double f) {
        if (g.containsKey(e)) {
            g.put(e, g.get(e) * f);
            h.put(e, h.get(e) + 1);
        } else {
            double min = 0;
            for (double d : g.values()) {
                if (min == 0 || d < min) {
                    min = d;
                }
            }
            if (min == 0) {
                min = 1; // the starting value
            }
            g.put(e, min);
            h.put(e, 1);
        }
    }

    private Lattice applyMove(Lattice l) {
        Lattice moved = new Lattice(l);
        // TODO: moves
        return moved;
    }

    @Override
    public Map<Double, Double> getDensity(int dimension, Polypeptide polypeptide) {
        int count = 0;
        double f = Math.E;
        g.clear();
        Lattice old = new Lattice(dimension);
        int size = polypeptide.size();
        for (int i=0; i<size; i++) {
            old.put(Point.point(i, 0, 0), polypeptide.get(i));
        }
        while (f > F_FINAL) {
            h.clear();
            while (!isSufficientlyFlat()) {
                Lattice trial = applyMove(old);
                if (accept(old, trial)) {
                    updateMaps(trial.getEnergy(), f);
                    old = trial;
                } else {
                    updateMaps(old.getEnergy(), f);
                }
                count++;
                if (count % 1000000 == 0) {
                    System.out.println((count / 1000000) + "M states counted");
                }
            }
            f = Math.sqrt(f);
        }
        return g;
    }
}
