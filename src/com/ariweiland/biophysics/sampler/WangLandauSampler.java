package com.ariweiland.biophysics.sampler;

import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.RandomUtils;
import com.ariweiland.biophysics.lattice.MovableLattice;
import com.ariweiland.biophysics.lattice.PullMove;
import com.ariweiland.biophysics.lattice.RebridgeMove;
import com.ariweiland.biophysics.peptide.Polypeptide;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Ari Weiland
 */
public class WangLandauSampler extends Sampler {

    public static final int MIN_SAMPLES = 20;
    public static final double FLATNESS = 0.8; // must be between 0 and 1 exclusive
    public static final double F_FINAL = 0.00000001; // 10^-8
    public static final double MOVE_RATIO = 0.2;

    private final Map<Double, Double> g = new HashMap<>();
    private final Map<Double, Integer> h = new HashMap<>();

    private boolean isSufficientlyFlat() {
        if (g.isEmpty() || h.size() != g.size()) {
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
                min = Double.MIN_VALUE; // the starting value
            }
            g.put(energy, min);
        }
        return g.get(energy);
    }

    @Override
    public Map<Double, Double> getDensity(int dimension, Polypeptide polypeptide) {
        int count = 0;
        int pullCount = 0;
        int rebridgeCount = 0;
        double f = Math.E;
        g.clear();
        int size = polypeptide.size();
        MovableLattice a = new MovableLattice(dimension, size);
        for (int i=0; i<size; i++) {
            a.put(new Point(i, 0, 0), polypeptide.get(i));
        }
        updateMaps(a.getEnergy(), f);
        while (Math.log(f) > F_FINAL) {
            h.clear();
            while (!isSufficientlyFlat()) {
                MovableLattice b = new MovableLattice(a);
                List<RebridgeMove> rebridgeMoves = a.getRebridgeMoves();
                boolean isPull = rebridgeMoves.isEmpty() || RandomUtils.tryChance(MOVE_RATIO);
                double threshold = 1;
                if (isPull) { // pull move
                    List<PullMove> aMoves = a.getPullMoves();
                    PullMove move = RandomUtils.selectRandom(aMoves);
                    b.pull(move);
                    List<PullMove> bMoves = b.getPullMoves();
                    threshold = ((double) aMoves.size()) / bMoves.size();
                } else {      // bond-rebridging move
                    RebridgeMove move = RandomUtils.selectRandom(rebridgeMoves);
                    b.rebridge(move);
                }
                threshold *= g(a.getEnergy()) / g(b.getEnergy());
                // potentially slightly faster because randoms do not always need to be generated
                if (RandomUtils.tryChance(threshold)) {
                    updateMaps(b.getEnergy(), f);
                    a = b;
                    if (isPull) {
                        pullCount++;
                    } else {
                        rebridgeCount++;
                    }
                } else {
                    updateMaps(a.getEnergy(), f);
                }
                count++;
                if (count % 100000 == 0) {
                    System.out.println((count / 1000) + "K trials");
                    System.out.println("\t" + (pullCount) + " pull moves");
                    System.out.println("\t" + (rebridgeCount) + " rebridge moves");
                }
            }
            f = Math.sqrt(f);
        }
        System.out.println(count + " total trials");
        System.out.println(pullCount + " total pull moves");
        System.out.println(rebridgeCount + " total rebridge moves");
        return g;
    }
}
