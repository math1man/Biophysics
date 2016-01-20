package com.ariweiland.biophysics.density;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.lattice.Folding;
import com.ariweiland.biophysics.lattice.Lattice;
import com.ariweiland.biophysics.peptide.Polypeptide;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Ari Weiland
 */
public class BruteForceHistogram extends Histogram {

    @Override
    public Map<Double, Integer> count(int dimension, Polypeptide polypeptide) {
        int size = polypeptide.size();
        Map<Double, Integer> counter = new HashMap<>();
        int[] state = new int[size]; // we won't actually use the 0 index
        Arrays.fill(state, -1);
        Lattice first = new Lattice(dimension, size);
        first.put(Point.point(0, 0, 0), polypeptide.get(0));
        first.put(Point.point(1, 0, 0), polypeptide.get(1));
        Folding[] foldings = new Folding[size];
        foldings[1] = new Folding(first, Point.point(1, 0, 0), 1, 0);
        int foldIndex = 2; // the first residue is directionless, and the second one just determines symmetry
        int count = 0;
        while (foldIndex > 1) {
            // Increment the position of the queen in the currently specified row
            state[foldIndex]++;
            // If we go over, reset that row and decrement the foldIndex
            // so that it will increment the previous row position
            // If we are at the first row, placing the queen past
            // half way is symmetric, so we can ignore them
            if (state[foldIndex] >= 2 * dimension) {
                state[foldIndex] = -1;
                foldIndex--;
            } else {
                Point next = foldings[foldIndex - 1].lastPoint.getAdjacent(Direction.values()[state[foldIndex]]);
                Lattice lattice = new Lattice(foldings[foldIndex - 1].lattice);
                // Check that the generated state is valid
                if (!lattice.containsPoint(next)) {
                    lattice.put(next, polypeptide.get(foldIndex));
                    if (foldIndex < size - 1) {
                        // If valid and there are still more rows to position, go to the next row
                        foldings[foldIndex] = new Folding(lattice, next, foldIndex, 0);
                        foldIndex++;
                    } else {
                        // Otherwise, increment the counter, reset the current row, and go back
                        double energy = lattice.getEnergy();
                        if (!counter.containsKey(energy)) {
                            counter.put(energy, 0);
                        }
                        counter.put(energy, 1 + counter.get(energy));
                        state[foldIndex] = -1;
                        foldIndex--;
                        count++;
                        if (count % 1000000 == 0) {
                            System.out.println(count + " states counted");
                        }
                    }
                }
            }
        }
        System.out.println(count + " total states counted");
        return counter;
    }
}
