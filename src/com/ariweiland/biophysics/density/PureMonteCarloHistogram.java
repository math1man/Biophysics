package com.ariweiland.biophysics.density;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.lattice.Lattice;
import com.ariweiland.biophysics.peptide.Polypeptide;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Ari Weiland
 */
public class PureMonteCarloHistogram extends Histogram {

    private final int samples;

    public PureMonteCarloHistogram(int samples) {
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
                Point next;
                do {
                    Direction d = Direction.values()[((int) (Math.random() * 2 * dimension))];
                    next = last.getAdjacent(d);
                } while (lattice.containsPoint(next));
                lattice.put(next, polypeptide.get(j));
                last = next;
                isBoxedIn = true;
                for (Direction d : Direction.values(dimension)) {
                    if (!lattice.containsPoint(last.getAdjacent(d))) {
                        isBoxedIn = false;
                    }
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
                    System.out.println(count + " states counted");
                }
            }
        }
        System.out.println(count + " total states counted");
        return counter;
    }
}
