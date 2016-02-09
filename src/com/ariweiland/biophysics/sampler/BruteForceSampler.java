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
public class BruteForceSampler extends Sampler {

    @Override
    public Map<Double, Double> getDensity(int dimension, Polypeptide polypeptide) {
        int size = polypeptide.size();
        Map<Double, Double> counter = new HashMap<>();
        int[] state = new int[size]; // we won't actually use the 0 index
        Arrays.fill(state, -1);
        BacktrackLattice lattice = new BacktrackLattice(dimension, size);
        lattice.put(new Point(0, 0, 0), polypeptide.get(0));
        lattice.put(new Point(1, 0, 0), polypeptide.get(1));
        long count = 0;
        while (lattice.size() > 1) {
            // Change the direction of the currently specified residue
            int index = lattice.size();
            state[index]++;
            // If we go over, reset that direction and remove the last residue
            if (state[index] >= 2 * dimension) {
                state[index] = -1;
                lattice.removeLast();
            } else {
                Point next = lattice.getLastPoint().getAdjacent(Direction.values()[state[index]]);
                // Check that the generated state is valid
                if (!lattice.containsPoint(next)) {
                    lattice.put(next, polypeptide.get(index));
                    if (lattice.size() == size) {
                        // Otherwise, increment the counter
                        double energy = lattice.getEnergy();
                        if (!counter.containsKey(energy)) {
                            counter.put(energy, 0.0);
                        }
                        counter.put(energy, 1 + counter.get(energy));
                        // and remove the last residue
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
