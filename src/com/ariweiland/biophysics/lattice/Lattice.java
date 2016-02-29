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
    private final boolean hasSurface;
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
        this.hasSurface = (surface != null);
        this.surface = surface;
        this.lattice = new HashMap<>();
    }

    public Lattice(int dimension, int initialCapacity) {
        this(dimension, initialCapacity, null);
    }

    public Lattice(int dimension, int initialCapacity, Residue surface) {
        if (dimension < 2 || dimension > 3) {
            throw new IllegalArgumentException("Dimension of less than 2 or more than 3 does not make sense");
        }
        this.dimension = dimension;
        this.hasSurface = (surface != null);
        this.surface = surface;
        this.lattice = new HashMap<>(initialCapacity);
    }

    public Lattice(Lattice lattice) {
        this.dimension = lattice.dimension;
        this.hasSurface = lattice.hasSurface;
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
        return hasSurface;
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
        return lattice.containsKey(point) || (hasSurface && point.y == 0);
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

    /**
     * Draws an ASCII visualization of the peptides in the lattice to the console.
     * Also returns the drawing as a list of strings. Currently only draws 2D lattices.
     */
    public List<String> visualize() {
        int plusXBound = 0;
        int minusXBound = 0;
        int plusYBound = 0;
        int minusYBound = 0;
        int plusZBound = 0;
        int minusZBound = 0;
        boolean isFirst = true;
        for (Point point : points()) {
            if (isFirst) {
                plusXBound  = point.x;
                minusXBound = point.x;
                plusYBound  = point.y;
                minusYBound = point.y;
                plusZBound  = point.z;
                minusZBound = point.z;
                isFirst = false;
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
        }
        List<String> lines = new ArrayList<>();
        for (int k=minusZBound; k<=plusZBound; k++) {
            for (int i=plusYBound; i>=minusYBound; i--) {
                StringBuilder latticeString = new StringBuilder();
                StringBuilder connectionsString = new StringBuilder();
                for (int j=minusXBound; j<=plusXBound; j++) {
                    Peptide p = get(new Point(j, i, k));
                    if (p != null) {
                        int index = p.index;
                        String residue = p.residue.toString();
                        Point up = new Point(j, i, k + 1);
                        if (contains(up) && (get(up).index == index + 1 || get(up).index == index - 1)) {
                            residue = residue.replace('(', '{').replace(')', '}');
                        }
                        latticeString.append(residue);
                        Point right = new Point(j + 1, i, k);
                        if (contains(right) && (get(right).index == index + 1 || get(right).index == index - 1)) {
                            latticeString.append("-");
                        } else {
                            latticeString.append(" ");
                        }
                        Point below = new Point(j, i - 1, k);
                        if (contains(below) && (get(below).index == index + 1 || get(below).index == index - 1)) {
                            connectionsString.append(" |  ");
                        } else {
                            connectionsString.append("    ");
                        }
                    } else {
                        latticeString.append("    ");
                        connectionsString.append("    ");
                    }
                }
                lines.add(latticeString.toString());
                lines.add(connectionsString.toString());
                System.out.println(latticeString);
                System.out.println(connectionsString);
            }
            if (hasSurface()) {
                StringBuilder surface = new StringBuilder();
                StringBuilder base = new StringBuilder();
                for (int j=minusXBound; j<=plusXBound; j++) {
                    surface.append(getSurface()).append(" ");
                    base.append("-+--");
                }
                lines.add(surface.toString());
                lines.add(base.toString());
                System.out.println(surface);
                System.out.println(base);
            } else {
                StringBuilder edge = new StringBuilder();
                for (int j=minusXBound; j<=plusXBound; j++) {
                    edge.append("====");
                }
                lines.add(edge.toString());
                System.out.println(edge);
            }
        }
        return lines;
    }
}
