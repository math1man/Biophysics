package com.ariweiland.biophysics.lattice;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.peptide.Residue;

import java.util.*;

/**
 * This is the most basic form of a lattice. It keeps track of only the most essential
 * features, such as dimesion, surface, and energy. It also does not do any error checking
 * except that the dimension must be initialized as either 2 or 3. The contains, get, and
 * put methods do not have sanity checks.
 *
 * @author Ari Weiland
 */
public class Lattice {

    private final int dimension;
    private final Residue surface;
    protected final Map<Point, Peptide> lattice;
    protected double energy = 0;

    public Lattice(int dimension) {
        this(dimension, null);
    }

    public Lattice(int dimension, Residue surface) {
        if (dimension < 2 || dimension > 3) {
            throw new IllegalArgumentException("Dimension of less than 2 or more than 3 does not make sense");
        }
        this.dimension = dimension;
        this.lattice = new HashMap<>();
        this.surface = surface;
    }

    public Lattice(int dimension, int initialCapacity) {
        this(dimension, initialCapacity, null);
    }

    public Lattice(int dimension, int initialCapacity, Residue surface) {
        if (dimension < 2 || dimension > 3) {
            throw new IllegalArgumentException("Dimension of less than 2 or more than 3 does not make sense");
        }
        this.dimension = dimension;
        this.lattice = new HashMap<>(initialCapacity);
        this.surface = surface;
    }

    public Lattice(Lattice lattice) {
        this.dimension = lattice.dimension;
        this.surface = lattice.surface;
        this.lattice = new HashMap<>(lattice.lattice);
        this.energy = lattice.energy;
    }

    public int getDimension() {
        return dimension;
    }

    /**
     * Returns whether or not this lattice has a surface
     * @return
     */
    public boolean hasSurface() {
        return surface != null;
    }

    /**
     * Returns the surface of this lattice, or null if it does not have a surface
     * @return
     */
    public Residue getSurface() {
        return surface;
    }

    /**
     * Returns the number of peptides in the lattice
     * @return
     */
    public int size() {
        return lattice.size();
    }

    /**
     * Returns true if the lattice is empty
     * @return
     */
    public boolean isEmpty() {
        return lattice.isEmpty();
    }

    /**
     * Returns true if the specified point is occupied
     * @param point
     * @return
     */
    public boolean contains(Point point) {
        return lattice.containsKey(point) || (surface != null && point.y == 0);
    }

    /**
     * Returns the peptide at the specified point, or null.
     *
     * If there is a a surface and point.y == 0, returns a peptide whose residue is the
     * surface residue, and whose index is -2 so as to not interact with other peptides.
     *
     * @param point
     * @return
     */
    public Peptide get(Point point) {
        if (hasSurface() && point.y == 0) {
            return new Peptide(-2, getSurface());
        } else {
            return lattice.get(point);
        }
    }

    /**
     * Places the peptide in the lattice, and updates the energy appropriately.
     * If you place a peptide at a previously occupied point, it will NOT throw
     * an exception, and the lattice will be invalid. For a lattice with sanity
     * checking, use CheckedLattice.
     *
     * @param point
     * @param peptide
     */
    public void put(Point point, Peptide peptide) {
        for (Direction d : Direction.values(getDimension())) {
            Peptide adj = get(point.getAdjacent(d));
            if (adj != null) {
                // if they are not adjoining peptides
                if (adj.index != peptide.index + 1 && adj.index != peptide.index - 1) {
                    energy += peptide.interaction(adj);
                }
                energy -= adj.interaction(Residue.H2O);
            } else {
                energy += peptide.interaction(Residue.H2O);
            }
        }
        lattice.put(point, peptide);
    }

    /**
     * Puts the contents of another lattice map into this lattice
     * @param lattice
     */
    public void putAll(Map<Point, Peptide> lattice) {
        for (Map.Entry<Point, Peptide> e : lattice.entrySet()) {
            put(e.getKey(), e.getValue());
        }
    }

    /**
     * Clears the lattice and sets the energy to 0
     */
    public void clear() {
        lattice.clear();
        energy = 0;
    }

    /**
     * Returns a set of all occupied points in the lattice
     * @return
     */
    public Set<Point> points() {
        return lattice.keySet();
    }

    /**
     * Returns the lattice energy
     * @return
     */
    public double getEnergy() {
        return Math.round(energy * 100) / 100.0;
    }

}
