package com.ariweiland.biophysics.lattice;

import com.ariweiland.biophysics.Point;

/**
 * This class is a simple wrapper class for a partial folding of a polypeptide.
 * It contains a Lattice, the last point in the lattice a peptide was added, the
 * index of the last peptide to be added, and the energy bound of the lattice.
 *
 * @author Ari Weiland
 */
public class Folding implements Comparable<Folding> {

    public final Lattice lattice;
    public final Point lastPoint;
    public final int index;
    public final double energyBound;

    public Folding(Lattice lattice, int lastX, int lastY, int index, double energyBound) {
        this(lattice, new Point(lastX, lastY), index, energyBound);
    }

    public Folding(Lattice lattice, Point lastPoint, int index, double energyBound) {
        this.lattice = lattice;
        this.lastPoint = lastPoint;
        this.index = index;
        this.energyBound = energyBound;
    }

    @Override
    public int compareTo(Folding o) {
        int compare = Double.compare(energyBound, o.energyBound);
        if (compare == 0) {
            compare = Integer.compare(lattice.getSurfaceSize(), o.lattice.getSurfaceSize());
        }
        return compare;
    }

    @Override
    public String toString() {
        return energyBound + "/" + lattice.getSurfaceSize();
    }
}
