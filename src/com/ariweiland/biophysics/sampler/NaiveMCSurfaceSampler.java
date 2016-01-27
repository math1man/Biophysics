package com.ariweiland.biophysics.sampler;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.lattice.SurfaceLattice;
import com.ariweiland.biophysics.peptide.Polypeptide;
import com.ariweiland.biophysics.peptide.Residue;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Ari Weiland
 */
public class NaiveMCSurfaceSampler extends Sampler {

    private final int samples;
    private final Residue surface;

    public NaiveMCSurfaceSampler(int samples, Residue surface) {
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
    public Map<Double, Double> getDensity(int dimension, Polypeptide polypeptide) {
        int size = polypeptide.size();
        Map<Double, Double> counter = new HashMap<>();
        int count = 0;
        for (int i=0; i<samples; i++) {
            SurfaceLattice lattice = new SurfaceLattice(dimension, surface, size);
            // start out at a random y value between 1 and size, inclusive
            Point last = Point.point(0, (int) (Math.random() * size) + 1, 0);
            lattice.put(last, polypeptide.get(0));
            boolean isBoxedIn = false;
            for (int j=1; j<size && !isBoxedIn; j++) {
                List<Direction> opens = new ArrayList<>();
                for (Direction d : Direction.values(dimension)) {
                    if (!lattice.containsPoint(last.getAdjacent(d)) && last.getAdjacent(d).y <= size) {
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
