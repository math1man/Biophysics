package com.ariweiland.biophysics.sampler;

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
public class BruteForceSampler extends Sampler {

    @Override
    public Map<Double, Double> getDensity(int dimension, Polypeptide polypeptide) {
        int size = polypeptide.size();
        Map<Double, Double> counter = new HashMap<>();
        int[] state = new int[size]; // we won't actually use the 0 index
        Arrays.fill(state, -1);
        Lattice first = new Lattice(dimension, size);
        first.put(new Point(0, 0, 0), polypeptide.get(0));
        first.put(new Point(1, 0, 0), polypeptide.get(1));
        Folding[] foldings = new Folding[size];
        foldings[1] = new Folding(first, new Point(1, 0, 0), 1, 0);
        int foldIndex = 2; // the first residue is directionless, and the second one just determines symmetry
        long count = 0;
        while (foldIndex > 1) {
            // Change the direction of the currently specified residue
            state[foldIndex]++;
            // If we go over, reset that direction and decrement the foldIndex
            // so that it will modify the previous residue direction
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
                        // If valid and there are still more residues to position, go to the next residue
                        foldings[foldIndex] = new Folding(lattice, next, foldIndex, 0);
                        foldIndex++;
                    } else {
                        // Otherwise, increment the counter, reset the current residues, and go back
                        double energy = lattice.getEnergy();
                        if (!counter.containsKey(energy)) {
                            counter.put(energy, 0.0);
                        }
                        counter.put(energy, 1 + counter.get(energy));
                        state[foldIndex] = -1;
                        foldIndex--;
                        count++;
                        if (count % 1000000 == 0) {
                            System.out.println((count / 1000000) + "M states counted");
                        }
                    }
                }
            }
        }
        System.out.println(count + " total states counted");
        return counter;
    }
}
