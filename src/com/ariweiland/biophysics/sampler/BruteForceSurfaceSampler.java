package com.ariweiland.biophysics.sampler;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.lattice.BacktrackLattice;
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
    private boolean running;

    public BruteForceSurfaceSampler(Residue surface) {
        this.surface = surface;
    }


    public Residue getSurface() {
        return surface;
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
        long count = 0;
        int maxY = size + 1;
        for (int y = 1; y < maxY && running; y++) {
            int[] state = new int[size]; // we won't actually use the 0 index
            Arrays.fill(state, -1);
            BacktrackLattice lattice = new BacktrackLattice(dimension, size, surface);
            lattice.put(new Point(0, y, 0), polypeptide.get(0));
            while (!lattice.isEmpty() && running) {
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
                    if (!lattice.contains(next) && next.y < maxY) { // TODO: consider next.y <= maxY
                        lattice.put(next, polypeptide.get(index));
                        if (lattice.size() == size) {
                            // Otherwise, increment the counter, reset the current row, and go back
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
        }
        System.out.println(count + " total states counted");
        return counter;
    }
}
