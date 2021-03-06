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
public class NaiveSampler extends Sampler {

    private int samples;

    private boolean running;
    public NaiveSampler() {
        this(10000000);
    }

    public NaiveSampler(int samples) {
        this.samples = samples;
    }

    public int getSamples() {
        return samples;
    }

    public void setSamples(int samples) {
        this.samples = samples;
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
                    Point next = last.getAdjacent(RandomUtils.selectRandom(opens));
                    lattice.put(next, polypeptide.get(j));
                    last = next;
                }
            }
            if (!isBoxedIn) {
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
