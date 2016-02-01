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
public class AriMCSampler extends Sampler {

    private final int samples;
    private final double minChance = 1.0;

    public AriMCSampler(int samples) {
        this.samples = samples;
    }

    public int getSamples() {
        return samples;
    }

    /**
     * This method calculates the angle betwen the vector from a to b and a to (0,0,0).
     * If the angle is 0, there is a 50% chance of acceptance, and if it is 180, there
     * is an 100% chance of acceptance. Angles in between scale linearly.
     * @param a
     * @param b
     * @return
     */
    private boolean accept1(Point a, Point b) {
        int x = b.x - a.x;
        int y = b.y - a.y;
        int z = b.z - a.z;
        // theta has range 0 to Pi
        double theta = Math.acos(Math.sqrt(x * b.x + y * b.y + z * b.z)
                / Math.sqrt(b.x * b.x + b.y * b.y + b.z * b.z)
                / Math.sqrt(x * x + y * y + z * z));
        // adjust range from 0.5 to 1
        double test = (1 - minChance) * theta / Math.PI + minChance;
        return RandomUtils.tryChance(test);
    }

    /**
     * This method calculates the angle betwen the vector from a to b and a to (0,0,0).
     * If the angle is 0, there is a 50% chance of acceptance, and if it is 180, there
     * is an 100% chance of acceptance. Angles in between scale linearly.
     * @param a
     * @param b
     * @return
     */
    private boolean accept2(Point a, Point b) {
        int x = b.x - a.x;
        int y = b.y - a.y;
        int z = b.z - a.z;
        // cos has range 1 to -1
        double cos = Math.sqrt(x * b.x + y * b.y + z * b.z)
                / Math.sqrt(b.x * b.x + b.y * b.y + b.z * b.z)
                / Math.sqrt(x * x + y * y + z * z);
        double test = 0.5 + minChance / 2 - cos / (1 - minChance);
        return RandomUtils.tryChance(test);
    }

    @Override
    public Map<Double, Double> getDensity(int dimension, Polypeptide polypeptide) {
        int size = polypeptide.size();
        Map<Double, Double> counter = new HashMap<>();
        Lattice base = new Lattice(dimension, size);
        base.put(new Point(0, 0, 0), polypeptide.get(0));
        base.put(new Point(1, 0, 0), polypeptide.get(1));
        int count = 0;
        for (int i=0; i<samples; i++) {
            Lattice lattice = new Lattice(base);
            Point last = new Point(1, 0, 0);
            boolean isBoxedIn = false;
            for (int j=2; j<size && !isBoxedIn; j++) {
                List<Direction> opens = new ArrayList<>();
                for (Direction d : Direction.values(dimension)) {
                    if (!lattice.containsPoint(last.getAdjacent(d))) {
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
                    } while (!accept1(last, next));
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
                if (count % 1000000 == 0) {
                    System.out.println((count / 1000000) + "M states counted");
                }
            }
        }
        System.out.println(count + " total states counted");
        return counter;
    }
}
