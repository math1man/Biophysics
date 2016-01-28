package com.ariweiland.biophysics.sampler;

import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.lattice.MovableLattice;
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

    private boolean accept(MovableLattice old, MovableLattice trial) {
        double threshold = g(old.getEnergy()) / g(trial.getEnergy());
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
        g.put(e, g(e) * f);
        if (h.containsKey(e)) {
            h.put(e, h.get(e) + 1);
        } else {
            h.put(e, 1);
        }
    }

    private double g(double energy) {
        if (!g.containsKey(energy)) {
            double min = 0;
            for (double d : g.values()) {
                if (min == 0 || d < min) {
                    min = d;
                }
            }
            if (min == 0) {
                min = 1; // the starting value
            }
            g.put(energy, min);
        }
        return g.get(energy);
    }

    private MovableLattice applyMove(MovableLattice l) {
        MovableLattice moved = new MovableLattice(l);
        boolean success = false;
        while (!success) {
            int i = randomInt(moved.size() - 1);
            if (Math.random() < 0.2) { // pull move
                success = moved.pull(i);
            } else {                   // bond-rebridging move
                success = moved.rebridge(i);
            }
        }
        return moved;
    }

    private static int randomInt(int max) {
        return (int) (Math.random() * max);
    }

    @Override
    public Map<Double, Double> getDensity(int dimension, Polypeptide polypeptide) {
        int count = 0;
        double f = Math.E;
        g.clear();
        int size = polypeptide.size();
        MovableLattice old = new MovableLattice(dimension, size);
        for (int i=0; i<size; i++) {
            old.put(new Point(i, 0, 0), polypeptide.get(i));
        }
        updateMaps(old.getEnergy(), f);
        while (f > F_FINAL) {
            h.clear();
            while (!isSufficientlyFlat()) {
                MovableLattice trial = applyMove(old);
                if (accept(old, trial)) {
                    updateMaps(trial.getEnergy(), f);
                    old = trial;
                } else {
                    updateMaps(old.getEnergy(), f);
                }
                count++;
                if (count % 1000000 == 0) {
                    System.out.println((count / 1000000) + "M trials");
                }
            }
            f = Math.sqrt(f);
        }
        System.out.println(count + " total trials");
        return g;
    }
}
