package com.ariweiland.biophysics.lattice;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.peptide.Residue;

/**
 * This class represents a folded polypeptide, laid out in 2D rectangular grid.
 * It also keeps track of various relevant properties, such as the bounding box/cube
 * and perimeter/surface area of the polypeptide in the lattice and the lattice energy.
 * @author Ari Weiland
 */
public class CheckedLattice extends Lattice {

    protected double energy = 0;
    private int surfaceSize = 0;
    private int plusXBound = 0;
    private int minusXBound = 0;
    private int plusYBound = 0;
    private int minusYBound = 0;
    private int plusZBound = 0;
    private int minusZBound = 0;

    public CheckedLattice(int dimension) {
        this(dimension, null);
    }

    public CheckedLattice(int dimension, Residue surface) {
        super(dimension, surface);
    }

    public CheckedLattice(int dimension, int initialCapacity) {
        this(dimension, initialCapacity, null);
    }

    public CheckedLattice(int dimension, int initialCapacity, Residue surface) {
        super(dimension, initialCapacity, surface);
    }

    public CheckedLattice(CheckedLattice lattice) {
        super(lattice);
        this.surfaceSize = lattice.surfaceSize;
        this.plusXBound  = lattice.plusXBound;
        this.minusXBound = lattice.minusXBound;
        this.plusYBound  = lattice.plusYBound;
        this.minusYBound = lattice.minusYBound;
        this.plusZBound  = lattice.plusZBound;
        this.minusZBound = lattice.minusZBound;
    }

    @Override
    public void put(Point point, Peptide peptide) {
        if (getDimension() == 2 && point.z != 0) {
            throw new IllegalArgumentException("2D points cannot have a z-component");
        }
        if (hasSurface() && point.y < 1) {
            throw new IllegalArgumentException("Cannot put a point on or below the surface (y <= 0)");
        }
        if (contains(point)) {
            throw new IllegalArgumentException("That point is already occupied");
        }
        if (isEmpty()) {
            plusXBound  = point.x;
            minusXBound = point.x;
            plusYBound  = point.y;
            minusYBound = point.y;
            plusZBound  = point.z;
            minusZBound = point.z;
        } else {
            if (point.x > plusXBound) {
                plusXBound = point.x;
            } else if (point.x < minusXBound) {
                minusXBound = point.x;
            }
            if (point.y > plusYBound) {
                plusYBound = point.y;
            } else if (point.y < minusYBound) {
                minusYBound = point.y;
            }
            if (point.z > plusZBound) {
                plusZBound = point.z;
            } else if (point.z < minusZBound) {
                minusZBound = point.z;
            }
        }
        for (Direction d : Direction.values(getDimension())) {
            Peptide adj = get(point.getAdjacent(d));
            if (adj == null) {
                surfaceSize += 1;
            } else if (adj.index >= 0) { // this is for use with surface lattice, so that it properly handles surface-perimeter
                surfaceSize -= 1;
            }
        }
        super.put(point, peptide);
    }

    @Override
    public void clear() {
        super.clear();
        surfaceSize = 0;
        plusXBound  = 0;
        minusXBound = 0;
        plusYBound  = 0;
        minusYBound = 0;
        plusZBound  = 0;
        minusZBound = 0;
    }

    /**
     * Returns the surface size of the peptides in the lattice
     * @return
     */
    public int getSurfaceSize() {
        return surfaceSize;
    }

    /**
     * Returns the perimeter of the smallest bounding box or surface area of
     * the smallest bounding cube that contains the peptides in the lattice
     * @return
     */
    public int boundingPerimeter() {
        int xRange = plusXBound - minusXBound + 1;
        int yRange = plusYBound - minusYBound + 1;
        if (getDimension() == 2) {
            return 2 * (xRange + yRange);
        } else {
            int zRange = plusZBound - minusZBound + 1;
            return 2 * (xRange * yRange + xRange * zRange + yRange * zRange);
        }
    }

}
