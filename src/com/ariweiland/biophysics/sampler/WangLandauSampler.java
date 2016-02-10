package com.ariweiland.biophysics.sampler;

import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.RandomUtils;
import com.ariweiland.biophysics.lattice.MovableLattice;
import com.ariweiland.biophysics.lattice.PullMove;
import com.ariweiland.biophysics.peptide.Polypeptide;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Ari Weiland
 */
public class WangLandauSampler extends Sampler {

    public static final MathContext MC = MathContext.DECIMAL32;
    public static final double F_FINAL = 0.00000001; // 10^-8

    private double minEnergy = 0;   // must be less than or equal to 0
    private double flatness = 0.8;  // must be between 0 and 1 exclusive
    private int moveCount = 1;      // must be positive
    private double moveRatio = 0.2; // must be between 0 and 1 exclusive

    private final Map<Double, BigDecimal> g = new HashMap<>();
    private final Map<Double, Integer> h = new HashMap<>();

    public WangLandauSampler() {}

    public WangLandauSampler(double minEnergy) {
        this.minEnergy = minEnergy;
    }

    public WangLandauSampler(double minEnergy, double flatness) {
        this.minEnergy = minEnergy;
        this.flatness = flatness;
    }

    public WangLandauSampler(double minEnergy, double flatness, int moveCount, double moveRatio) {
        this.minEnergy = minEnergy;
        this.flatness = flatness;
        this.moveCount = moveCount;
        this.moveRatio = moveRatio;
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

    public int getMoveCount() {
        return moveCount;
    }

    public void setMoveCount(int moveCount) {
        this.moveCount = moveCount;
    }

    public double getMoveRatio() {
        return moveRatio;
    }

    public void setMoveRatio(double moveRatio) {
        this.moveRatio = moveRatio;
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

    private void updateMaps(double energy, BigDecimal f) {
        g.put(energy, g(energy).multiply(f, MC));
        if (h.containsKey(energy)) {
            h.put(energy, h.get(energy) + 1);
        } else {
            h.put(energy, 1);
        }
    }

    private BigDecimal g(double energy) {
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

    @Override
    public Map<Double, Double> getDensity(int dimension, Polypeptide polypeptide) {
        int count = 0;
        int pullCount = 0;
        int rebridgeCount = 0;
        double f = Math.E;
        BigDecimal bdf = BigDecimal.valueOf(f);
        g.clear();
        if (minEnergy != 0) { // if there is a minimum energy specified, initialize g
            for (double e = minEnergy; e <= 0; e++) {
                g.put(e, BigDecimal.ZERO);
            }
        }
        int size = polypeptide.size();
        MovableLattice old = new MovableLattice(dimension, size);
        for (int i=0; i<size; i++) {
            old.put(new Point(i, 0, 0), polypeptide.get(i));
        }
        updateMaps(old.getEnergy(), bdf);
        while (Math.log(f) > F_FINAL) {
            h.clear();
            while (!isSufficientlyFlat()) {
                MovableLattice trial = new MovableLattice(old);
                List<PullMove> pullMoves = old.getPullMoves();
                int nOld = pullMoves.size();
                int pulls = 0;
                int rebridges = 0;
                for (int i=0; i<moveCount; i++) {
                    if (RandomUtils.tryChance(moveRatio)) { // pull move
                        PullMove move = RandomUtils.selectRandom(pullMoves);
                        trial.pull(move);
                        pulls++;
                        pullMoves = trial.getPullMoves();
                    } else if (trial.rebridge(RandomUtils.randomInt(trial.size()))) {
                        rebridges++;
                        pullMoves = trial.getPullMoves();
                    } else {
                        i--;
                        // nothing changed, so don't need to update pullMoves
                    }
                }
                double threshold = g(old.getEnergy()).divide(g(trial.getEnergy()), MC)
                        .multiply(BigDecimal.valueOf(((double) nOld) / pullMoves.size()), MC).doubleValue();
                // potentially slightly faster because randoms do not always need to be generated
                if (RandomUtils.tryChance(threshold)) {
                    old = trial;
                    pullCount += pulls;
                    rebridgeCount += rebridges;
                }
                updateMaps(old.getEnergy(), bdf);
                count++;
                if (count % 100000 == 0) {
                    System.out.println((count / 1000) + "K trials");
                    System.out.println("\t" + (pullCount) + " pull moves");
                    System.out.println("\t" + (rebridgeCount) + " rebridge moves");
                }
            }
            f = Math.sqrt(f);
            bdf = BigDecimal.valueOf(f);
        }
        System.out.println(count + " total trials");
        System.out.println(pullCount + " total pull moves");
        System.out.println(rebridgeCount + " total rebridge moves");

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
}
