package com.ariweiland.biophysics.lattice;

/**
 * This class is a simple wrapper class for a partial folding of a polypeptide.
 * It contains a Lattice, the last point in the lattice a peptide was added, the
 * index of the last peptide to be added, and the energy bound of the lattice.
 *
 * @author Ari Weiland
 */
public class Folding implements Comparable<Folding> {

    public final Lattice2D lattice;
    public final Point2D lastPoint;
    public final int index;
    public final double energyBound;

    public Folding(Lattice2D lattice, int lastX, int lastY, int index, double energyBound) {
        this(lattice, new Point2D(lastX, lastY), index, energyBound);
    }

    public Folding(Lattice2D lattice, Point2D lastPoint, int index, double energyBound) {
        this.lattice = lattice;
        this.lastPoint = lastPoint;
        this.index = index;
        this.energyBound = energyBound;
    }

    @Override
    public int compareTo(Folding o) {
        int compare = Double.compare(energyBound, o.energyBound);
        if (compare == 0) {
            compare = Integer.compare(lattice.getPerimeter(), o.lattice.getPerimeter());
        }
        return compare;
    }

    @Override
    public String toString() {
        return energyBound + "/" + lattice.getPerimeter();
    }
}
