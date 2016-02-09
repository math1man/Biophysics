package com.ariweiland.biophysics.modeler;

import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.lattice.Folding;
import com.ariweiland.biophysics.lattice.BoundingLattice;
import com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.peptide.Polypeptide;
import com.ariweiland.biophysics.peptide.Residue;

import java.util.concurrent.PriorityBlockingQueue;

/**
 * @author Ari Weiland
 */
public abstract class SurfaceModeler extends ParallelModeler {

    private final Residue surface;

    protected SurfaceModeler(int dimension, Residue surface) {
        super(dimension);
        this.surface = surface;
    }

    public Residue getSurface() {
        return surface;
    }

    /**
     * This helper method determines the initial minimum energy of a polypeptide that
     * is about to begin the folding process.
     * @param polypeptide
     * @return
     */
    protected abstract double getInitialEnergyBound(Polypeptide polypeptide);

    /**
     * This helper method calculates the modification to the lower energy bound in the
     * initial seeding stage of filling the heap. It therefore assumes that no peptide
     * is adjacent to any other peptide in the polypeptide, though they may be adjacent
     * to the surface.
     * @param y
     * @param p
     * @return
     */
    protected abstract double getBoundAdjust(int y, Peptide p);

    @Override
    protected int getSeedCount(Polypeptide polypeptide) {
        return super.getSeedCount(polypeptide) * 10;
    }

    @Override
    protected PriorityBlockingQueue<Folding> initializeHeap(Polypeptide polypeptide) {
        PriorityBlockingQueue<Folding> initialHeap = new PriorityBlockingQueue<>();
        int size = polypeptide.size();
        // use this so that we don't bother with the peptide floating far away from the surface
        int maxY = getMaxY(polypeptide);
        // fill the queue initially.  this avoids symmetrical solutions
        for (int i = 1; i < maxY; i++) {
            for (int j = 1; j < maxY; j++) {
                BoundingLattice lattice = new BoundingLattice(getDimension(), surface, size);
                double bound = getInitialEnergyBound(polypeptide);
                int k;
                // add some number of residues between 0 and all of them in a vertical line, either rising or falling
                for (k = 0; k <= Math.abs(i - j) && k < size; k++) {
                    Peptide next = polypeptide.get(k);
                    int y;
                    if (i > j) {
                        y = i - k;
                    } else {
                        y = i + k;
                    }
                    lattice.put(new Point(0, y), next);
                    bound += getBoundAdjust(y, next);
                }
                int lastX = 0;
                // if there is at least one residue left, add it to the right of the last residue
                if (k < size) {
                    Peptide next = polypeptide.get(k);
                    lastX = 1;
                    lattice.put(new Point(lastX, j), next);
                    bound += getBoundAdjust(j, next);
                }
                // if all residues have been placed, replace the bound with the actual lattice energy
                if (k >= size - 1) {
                    bound = lattice.getEnergy();
                }
                // add the lattice to the heap as a Folding
                initialHeap.add(new Folding(lattice, new Point(lastX, j), k, bound));
            }
        }
        return initialHeap;
    }

    /**
     * This method calculates the maximum y-value a polypeptide should ever reach in surface modeling.
     * It is related to the perimeter bound.
     * @param polypeptide
     * @return
     */
    protected int getMaxY(Polypeptide polypeptide) {
        return getSurfaceBound(polypeptide) / 4 + 2;
    }
}
