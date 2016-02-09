package com.ariweiland.biophysics.sampler;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.lattice.Folding;
import com.ariweiland.biophysics.lattice.CheckedLattice;
import com.ariweiland.biophysics.peptide.Polypeptide;
import com.ariweiland.biophysics.peptide.Residue;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Ari Weiland
 */
public class BruteForceSurfaceSampler extends Sampler {

    private final Residue surface;

    public BruteForceSurfaceSampler(Residue surface) {
        this.surface = surface;
    }

    public Residue getSurface() {
        return surface;
    }

    @Override
    public Map<Double, Double> getDensity(int dimension, Polypeptide polypeptide) {
        int size = polypeptide.size();
        Map<Double, Double> counter = new HashMap<>();
        long count = 0;
        int maxY = size + 1;
        for (int y = 1; y < maxY; y++) {
            int[] state = new int[size]; // we won't actually use the 0 index
            Arrays.fill(state, -1);
            CheckedLattice first = new CheckedLattice(dimension, size, surface);
            first.put(new Point(0, y, 0), polypeptide.get(0));
            // TODO: update this shit
            Folding[] foldings = new Folding[size];
            foldings[0] = new Folding(first, new Point(0, y, 0), 0, 0);
            int foldIndex = 1; // the first residue is directionless
            while (foldIndex > 0) {
                // Change the direction of the currently specified residue
                state[foldIndex]++;
                // If we go over, reset that direction and decrement the foldIndex
                // so that it will modify the previous residue direction
                if (state[foldIndex] >= 2 * dimension) {
                    state[foldIndex] = -1;
                    foldIndex--;
                } else {
                    Point next = foldings[foldIndex - 1].lastPoint.getAdjacent(Direction.values()[state[foldIndex]]);
                    CheckedLattice lattice = new CheckedLattice(foldings[foldIndex - 1].lattice);
                    // Check that the generated state is valid
                    if (!lattice.contains(next) && next.y < maxY) {
                        lattice.put(next, polypeptide.get(foldIndex));
                        if (foldIndex < size - 1) {
                            // If valid and there are still more rows to position, go to the next row
                            foldings[foldIndex] = new Folding(lattice, next, foldIndex, 0);
                            foldIndex++;
                        } else {
                            // Otherwise, increment the counter, reset the current row, and go back
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
                }
            }
        }
        System.out.println(count + " total states counted");
        return counter;
    }
}
