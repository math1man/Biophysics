package com.ariweiland.biophysics.lattice;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.peptide.Peptide;

import java.util.ArrayList;
import java.util.List;

/**
 * The standard lattice does not allow modifications to the lattice after a Residue has been dropped,
 * aside from the clear method. This lattice will have methods for both pull moves and bond-rebridging
 * moves. This lattice also does support surface size or bounding perimeter. It also dynamically
 * calculates energy instead of maintaining it throughout.
 * @author Ari Weiland
 */
public class MovableLattice extends Lattice {

    private final List<Point> pointSequence;
    private boolean modified = false;

    public MovableLattice(int dimension) {
        super(dimension);
        pointSequence = new ArrayList<>();
    }

    public MovableLattice(int dimension, int initialCapacity) {
        super(dimension, initialCapacity);
        pointSequence = new ArrayList<>(initialCapacity);
    }

    public MovableLattice(MovableLattice lattice) {
        super(lattice);
        pointSequence = new ArrayList<>(lattice.pointSequence);
    }

    @Override
    public void put(Point point, Peptide peptide) {
        if (getDimension() == 2 && point.z != 0) {
            throw new IllegalArgumentException("2D points cannot have a z-component");
        }
        if (containsPoint(point)) {
            throw new IllegalArgumentException("That point is already occupied");
        }
        modified = true;
        lattice.put(point, peptide);
        pointSequence.add(point);
    }

    @Override
    public void clear() {
        super.clear();
        pointSequence.clear();
    }

    @Override
    public double getEnergy() {
        if (modified) { // lazy update
            energy = 0;
            for (Point p : pointSequence) {
                Peptide peptide = get(p);
                for (Direction d : Direction.values(getDimension())) {
                    Peptide adj = get(p.getAdjacent(d));
                    // don't count if lower index to avoid double counting
                    // don't count if the next index because that is not an interaction
                    if (adj == null || adj.index > peptide.index + 1) {
                        energy += peptide.interaction(adj);
                    }
                }
            }
        }
        return energy;
    }

    @Override
    public int getSurfaceSize() {
        throw new UnsupportedOperationException();
    }

    @Override
    public int boundingPerimeter() {
        throw new UnsupportedOperationException();
    }

    /**
     * Applies a pull move on the residue at index i to a random valid location.
     * Returns true if successful, or false if no valid location was found.
     * Throws an IllegalArgumentException if i specifies the end residue.
     * @param i
     * @return
     */
    public boolean pull(int i) {
        if (i == size() - 1) {
            throw new IllegalArgumentException("Cannot select the end residue");
        }
        modified = true;
        Point point = pointSequence.get(i);
        Point next = pointSequence.get(i + 1);
        List<Direction> options = new ArrayList<>();
        for (Direction d : Direction.values(getDimension())) {
            Point c = point.getAdjacent(d);
            Point l = next.getAdjacent(d);
            if (!containsPoint(l) // position L is open
                    && !next.getAdjacent(d.getReverse()).equals(point) // this is not collinear with point and next
                    && (!containsPoint(c) || get(c).index == i - 1)) { // point c is open or occupied by peptide i-1
                options.add(d);
            }
        }
        Direction d;
        if (options.isEmpty()) {
            return false;
        } else if (options.size() == 1) {
            d = options.get(0);
        } else {
            d = options.get((int) (Math.random() * options.size()));
        }
        Point c = point.getAdjacent(d);
        Point l = next.getAdjacent(d);
        lattice.put(l, lattice.remove(point)); // first move point to L
        pointSequence.add(i, l);
        if (!containsPoint(c) && i > 0) {
            int j = i - 1;
            next = l;           // the new location of point j + 1 (used to check adjacency)
            Point two = c;      // the former location of point j + 2 and the new location of point j
            Point one = point;  // the former location of point j + 1 (temporary variable)
            while (j >= 0) {
                point = pointSequence.get(j); // get the next point
                if (!point.isAdjacentTo(next)) {
                    lattice.put(two, lattice.remove(point)); // move point to two
                    pointSequence.add(j, point);
                    next = two;
                    two = one;
                    one = point;
                    j--;
                } else {
                    j = -1; // break
                }
            }
        }
        return true;
    }

    /**
     * Applies a bond-rebridging move...
     * @param index
     * @return
     */
    public boolean rebridge(int index) {
        modified = true;
        // TODO: bond rebridging moves
        return true;
    }
}
