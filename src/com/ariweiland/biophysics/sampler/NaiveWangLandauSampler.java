package com.ariweiland.biophysics.sampler;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.RandomUtils;
import com.ariweiland.biophysics.lattice.Lattice;
import com.ariweiland.biophysics.peptide.Polypeptide;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * This yields results pretty much identical to the basic naive Monte Carlo algorithm.
 * @author Ari Weiland
 */
public class NaiveWangLandauSampler extends Sampler {

    public static final double F_FINAL = 0.00000001; // 10^-8

    private final double flatness; // must be between 0 and 1 exclusive

    private final Map<Double, Double> g = new HashMap<>();
    private final Map<Double, Integer> h = new HashMap<>();

    public NaiveWangLandauSampler(double flatness) {
        this.flatness = flatness;
    }

    public double getFlatness() {
        return flatness;
    }

    private boolean isSufficientlyFlat() {
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
                min = Double.MIN_VALUE; // the starting value
            }
            g.put(energy, min);
        }
        return g.get(energy);
    }

    @Override
    public Map<Double, Double> getDensity(int dimension, Polypeptide polypeptide) {
        int size = polypeptide.size();
        double f = Math.E;
        g.clear();

        Lattice base = new Lattice(dimension, size);
        base.put(new Point(0, 0, 0), polypeptide.get(0));
        base.put(new Point(1, 0, 0), polypeptide.get(1));
        int count = 0;
        while (Math.log(f) > F_FINAL) {
            h.clear();
            Lattice old = null;
            while (!isSufficientlyFlat()) {
                Lattice trial = new Lattice(base);
                Point last = new Point(1, 0, 0);
                boolean isBoxedIn = false;
                for (int j=2; j<size && !isBoxedIn; j++) {
                    List<Direction> opens = new ArrayList<>();
                    for (Direction d : Direction.values(dimension)) {
                        if (!trial.containsPoint(last.getAdjacent(d))) {
                            opens.add(d);
                        }
                    }
                    if (opens.isEmpty()) {
                        isBoxedIn = true;
                    } else {
                        Direction d = RandomUtils.selectRandom(opens);
                        Point next = last.getAdjacent(d);
                        trial.put(next, polypeptide.get(j));
                        last = next;
                    }
                }
                if (trial.size() == size) {
                    if (old == null || RandomUtils.tryChance(g(old.getEnergy()) / g(trial.getEnergy()))) {
                        old = trial;
                    }
                    updateMaps(old.getEnergy(), f);
                    count++;
                    if (count % 1000000 == 0) {
                        System.out.println((count / 1000000) + "M trials");
                    }
                }
            }
            f = Math.sqrt(f);
        }
        System.out.println(count + " total states counted");
        return g;
    }
}
