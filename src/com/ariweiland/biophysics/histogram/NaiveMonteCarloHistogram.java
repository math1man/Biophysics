package com.ariweiland.biophysics.histogram;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.lattice.Lattice;
import com.ariweiland.biophysics.peptide.Polypeptide;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Ari Weiland
 */
public class NaiveMonteCarloHistogram extends Histogram {

    private final int samples;

    public NaiveMonteCarloHistogram(int samples) {
        this.samples = samples;
    }

    public int getSamples() {
        return samples;
    }

    @Override
    public Map<Double, Integer> count(int dimension, Polypeptide polypeptide) {
        int size = polypeptide.size();
        Map<Double, Integer> counter = new HashMap<>();
        Lattice base = new Lattice(dimension, size);
        base.put(Point.point(0, 0, 0), polypeptide.get(0));
        base.put(Point.point(1, 0, 0), polypeptide.get(1));
        int count = 0;
        for (int i=0; i<samples; i++) {
            Lattice lattice = new Lattice(base);
            Point last = Point.point(1, 0, 0);
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
                    Direction d = opens.get((int) (Math.random() * opens.size()));
                    Point next = last.getAdjacent(d);
                    lattice.put(next, polypeptide.get(j));
                    last = next;
                }
            }
            if (lattice.size() == size) {
                double energy = lattice.getEnergy();
                if (!counter.containsKey(energy)) {
                    counter.put(energy, 0);
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
