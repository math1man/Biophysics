package com.ariweiland.biophysics.histogram;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.lattice.SurfaceLattice;
import com.ariweiland.biophysics.peptide.Polypeptide;
import com.ariweiland.biophysics.peptide.Residue;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Ari Weiland
 */
public class PureMCSurfaceHistogram extends Histogram {

    private final int samples;
    private final Residue surface;

    public PureMCSurfaceHistogram(int samples, Residue surface) {
        this.samples = samples;
        this.surface = surface;
    }

    public int getSamples() {
        return samples;
    }

    public Residue getSurface() {
        return surface;
    }

    @Override
    public Map<Double, Integer> count(int dimension, Polypeptide polypeptide) {
        int size = polypeptide.size();
        Map<Double, Integer> counter = new HashMap<>();
        int count = 0;
        for (int i=0; i<samples; i++) {
            SurfaceLattice lattice = new SurfaceLattice(dimension, surface, size);
            // start out at a random y value between 1 and size, inclusive
            Point last = Point.point(0, (int) (Math.random() * size) + 1, 0);
            lattice.put(last, polypeptide.get(0));
            boolean isBoxedIn = false;
            for (int j=1; j<size && !isBoxedIn; j++) {
                Point next;
                do {
                    Direction d = Direction.values()[((int) (Math.random() * 2 * dimension))];
                    next = last.getAdjacent(d);
                } while (lattice.containsPoint(next) || next.y > size);
                lattice.put(next, polypeptide.get(j));
                last = next;
                isBoxedIn = true;
                for (Direction d : Direction.values(dimension)) {
                    if (!lattice.containsPoint(last.getAdjacent(d)) && last.getAdjacent(d).y <= size) {
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
                    System.out.println((count / 1000000) + "M states counted");
                }
            }
        }
        System.out.println(count + " total states counted");
        return counter;
    }
}
