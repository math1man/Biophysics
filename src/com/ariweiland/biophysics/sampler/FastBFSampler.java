package com.ariweiland.biophysics.sampler;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.lattice.BacktrackLattice;
import com.ariweiland.biophysics.peptide.Polypeptide;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Ari Weiland
 */
public class FastBFSampler extends Sampler {

    @Override
    public Map<Double, Double> getDensity(int dimension, Polypeptide polypeptide) {
        int size = polypeptide.size();
        Map<Double, Double> counter = new HashMap<>();
        int[] state = new int[size]; // we won't actually use the 0 index
        Arrays.fill(state, -1);
        BacktrackLattice lattice = new BacktrackLattice(dimension, size);
        lattice.put(new Point(0, 0, 0), polypeptide.get(0));
        lattice.put(new Point(1, 0, 0), polypeptide.get(1));
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
                lattice.removeLast();
            } else {
                Point next = lattice.getLastPoint().getAdjacent(Direction.values()[state[foldIndex]]);
                // Check that the generated state is valid
                if (!lattice.containsPoint(next)) {
                    lattice.put(next, polypeptide.get(foldIndex));
                    if (foldIndex < size - 1) {
                        // If valid and there are still more residues to position, go to the next residue
                        foldIndex++;
                    } else {
                        // Otherwise, increment the counter, reset the current residues, and go back
                        double energy = lattice.getEnergy();
                        if (!counter.containsKey(energy)) {
                            counter.put(energy, 0.0);
                        }
                        counter.put(energy, 1 + counter.get(energy));
                        lattice.removeLast();
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
