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
 * @author Ari Weiland
 */
public class WeilandNaiveSampler extends Sampler {

    private final int samples;
    private final double minChance;
    private boolean running;

    public WeilandNaiveSampler(int samples, double minChance) {
        this.samples = samples;
        this.minChance = minChance;
    }

    public int getSamples() {
        return samples;
    }

    /**
     * This method calculates the angle between the vector from a to b and a to (0,0,0).
     * If the angle is 0, there is a 50% chance of acceptance, and if it is 180, there
     * is an 100% chance of acceptance. Angles in between scale linearly.
     * @param a
     * @param b
     * @return
     */
    private boolean accept(Point a, Point b) {

        double a0 = a.x * a.x + a.y * a.y + a.z * a.z;
        double b0 = b.x * b.x + b.y * b.y + b.z * b.z;
        double ab = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z);

        double cos = (a0 + ab - b0) / (2 * Math.sqrt(a0) * Math.sqrt(ab));
        // theta has range 0 to Pi
        double theta = Math.acos(cos);
        // adjust range from 0.5 to 1
        double test = (1 - minChance) * theta / Math.PI + minChance;
        return RandomUtils.tryChance(test);
    }

    @Override
    public void terminate() {
        running = false;
    }

    @Override
    public Map<Double, Double> getDensity(int dimension, Polypeptide polypeptide) {
        running = true;
        int size = polypeptide.size();
        Map<Double, Double> counter = new HashMap<>();
        Lattice base = new Lattice(dimension, size);
        base.put(new Point(0, 0, 0), polypeptide.get(0));
        base.put(new Point(1, 0, 0), polypeptide.get(1));
        int count = 0;
        for (int i=0; i<samples && running; i++) {
            Lattice lattice = new Lattice(base);
            Point last = new Point(1, 0, 0);
            boolean isBoxedIn = false;
            for (int j=2; j<size && !isBoxedIn; j++) {
                List<Direction> opens = new ArrayList<>();
                for (Direction d : Direction.values(dimension)) {
                    if (!lattice.contains(last.getAdjacent(d))) {
                        opens.add(d);
                    }
                }
                if (opens.isEmpty()) {
                    isBoxedIn = true;
                } else {
                    Point next;
                    do {
                        Direction d = RandomUtils.selectRandom(opens);
                        next = last.getAdjacent(d);
                    } while (!accept(last, next));
                    lattice.put(next, polypeptide.get(j));
                    last = next;
                }
            }
            if (lattice.size() == size) {
                double energy = lattice.getEnergy();
                if (!counter.containsKey(energy)) {
                    counter.put(energy, 0.0);
                }
                counter.put(energy, 1 + counter.get(energy));
                count++;
                if (count % 100000 == 0) {
                    System.out.println((count / 1000) + "k states counted");
                }
            }
        }
        System.out.println(count + " total states counted");
        return counter;
    }
}
