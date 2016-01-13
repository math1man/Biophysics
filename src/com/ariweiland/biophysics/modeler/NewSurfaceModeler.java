package com.ariweiland.biophysics.modeler;

import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.lattice.Folding;
import com.ariweiland.biophysics.lattice.Lattice;
import com.ariweiland.biophysics.lattice.SurfaceLattice;
import com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.peptide.Polypeptide;
import com.ariweiland.biophysics.peptide.Residue;

import java.util.Queue;
import java.util.concurrent.PriorityBlockingQueue;

/**
 * @author Ari Weiland
 */
public class NewSurfaceModeler extends Modeler {

    private final Residue surface;

    public NewSurfaceModeler(int seedCount, Residue surface) {
        super(seedCount);
        this.surface = surface;
    }

    public Residue getSurface() {
        return surface;
    }

    /**
     * This helper method returns the min interaction of the surface,
     * adjusted such that the water-surface interaction is treated as zero.
     * @return
     */
    private double getAdjustedSurfaceMinInteraction() {
        return surface.minInteraction() - surface.interaction(Residue.H2O);
    }

    /**
     * This is a helper method that calculates the modification to the lower energy bound
     * in the initial seeding stage of filling the heap. It therefore assumes that no
     * peptide is adjacent to any other peptide in the polypeptide, though they may be
     * adjacent to the surface.
     *
     * First, it adds at least one of the peptide's favorable water interaction, and
     * subtracts the peptide's min interactions. Then, if it is adjacent to the surface,
     * it adds the interaction with the surface and subtracts any favorable water-surface
     * interaction. Otherwise, it adds another of the peptide's favorable water interactions,
     * and adds any unfavorable water-surface interaction because this peptide cannot block
     * one of those.
     *
     * @param y
     * @param p
     * @return
     */
    private double getBoundAdjust(int y, Peptide p) {
        double boundAdjust = getFavorableWaterInteraction(p) - 2 * p.minInteraction();
        if (y == 1) {
            boundAdjust += p.interaction(surface) - getAdjustedSurfaceMinInteraction();
        } else {
            boundAdjust += getFavorableWaterInteraction(p);
        }
        return boundAdjust;
    }

    public Lattice fold(Polypeptide polypeptide) {
        PriorityBlockingQueue<Folding> initialHeap = new PriorityBlockingQueue<>();
        int size = polypeptide.size();
        // use this so that we don't bother with the peptide floating far away from the surface
        int maxY = getMaxY(size);
        // fill the queue initially.  this avoids symmetrical solutions
        for (int i=1; i<maxY; i++) {
            for (int j=1; j<maxY; j++) {
                SurfaceLattice lattice = new SurfaceLattice(surface);
                double lowerBound = polypeptide.getMinEnergy()
                        + getFavorableWaterInteraction(polypeptide.get(0))
                        + size * getAdjustedSurfaceMinInteraction();
                int k;
                // add some number of residues between 0 and all of them in a vertical line, either rising or falling
                for (k=0; k<=Math.abs(i-j) && k<size; k++) {
                    Peptide next = polypeptide.get(k);
                    int y;
                    if (i > j) {
                        y = i - k;
                    } else {
                        y = i + k;
                    }
                    lattice.put(0, y, next);
                    lowerBound += getBoundAdjust(y, next);
                }
                int lastX = 0;
                // if there is at least one residue left, add it to the right of the last residue
                if (k < size) {
                    Peptide next = polypeptide.get(k);
                    lastX = 1;
                    lattice.put(lastX, j, next);
                    lowerBound += getBoundAdjust(j, next);
                }
                // if all residues have been placed, replace the bound with the actual lattice energy
                if (k >= size - 1) {
                    lowerBound = lattice.getEnergy();
                }
                // add the lattice to the heap as a Folding
                initialHeap.add(new Folding(lattice, lastX, j, k, lowerBound));
            }
        }

        // iterate a few times to make the initial heap bigger
        int count = 0;
        while (count < getSeedCount()) {
            Folding solution = iterate(polypeptide, initialHeap);
            if (solution != null) {
                return solution.lattice;
            }
            count++;
        }

        int processors = Runtime.getRuntime().availableProcessors();
        return parallelize(polypeptide, initialHeap, processors, MAX_HEAP_SIZE / processors - 1);
    }

    @Override
    public Folding iterate(Polypeptide polypeptide, Queue<Folding> queue) {
        int size = polypeptide.size();
        Folding folding = queue.poll();
        int nextIndex = folding.index + 1;
        if (nextIndex < size) {
            Peptide p = polypeptide.get(nextIndex);
            // try to add the peptide in every direction
            for (Point.Direction d : Point.Direction.values()) {
                Point next = folding.lastPoint.getAdjacent(d);
                if (!folding.lattice.containsPoint(next) && next.y < getMaxY(size)) {
                    SurfaceLattice l = new SurfaceLattice((SurfaceLattice) folding.lattice);
                    l.put(next, p);
                    // set the bound from the previous bound, minus the min interactions for this peptide,
                    // minus one favorable water interaction which
                    double bound = folding.energyBound - 2 * p.minInteraction() - getFavorableWaterInteraction(p);
                    if (nextIndex < size - 1) {
                        for (Point.Direction d1 : Point.Direction.values()) {
                            if (d1 != d.getReverse()) {
                                if (l.containsPoint(next.getAdjacent(d1))) {
                                    Peptide adjacent = l.get(next.getAdjacent(d1));
                                    bound += p.interaction(adjacent) - getFavorableWaterInteraction(adjacent);
                                } else {
                                    bound += getFavorableWaterInteraction(p);
                                }
                            }
                        }
                        if (next.y > 1) {
                            bound -= getAdjustedSurfaceMinInteraction();
                        }
                    } else {
                        bound = l.getEnergy();
                    }
                    queue.add(new Folding(l, next, nextIndex, bound));
                }
            }
            return null;
        } else {
            return folding;
        }
    }

    public static int getMaxY(int n) {
        return getPerimBound(n) / 4 + 2;
    }


}
