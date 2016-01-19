package com.ariweiland.biophysics.lattice;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.peptide.Residue;

import java.util.*;

/**
 * This class represents a folded polypeptide, laid out in 2D rectangular grid.
 * It also keeps track of various relevant properties, such as the bounding box/cube
 * and perimeter/surface area of the polypeptide in the lattice and the lattice energy.
 * @author Ari Weiland
 */
public class Lattice {

    private final int dimension;
    private final Map<Point, Peptide> lattice;
    private double energy = 0;
    private int surfaceSize = 0;
    protected int plusXBound = 0;
    protected int minusXBound = 0;
    protected int plusYBound = 0;
    protected int minusYBound = 0;
    protected int plusZBound = 0;
    protected int minusZBound = 0;

    public Lattice(int dimension) {
        if (dimension < 2 || dimension > 3) {
            throw new IllegalArgumentException("Dimension of less than 2 or more than 3 does not make sense");
        }
        this.dimension = dimension;
        this.lattice = new HashMap<>();
    }

    public Lattice(int dimension, int initialCapacity) {
        if (dimension < 2 || dimension > 3) {
            throw new IllegalArgumentException("Dimension of less than 2 or more than 3 does not make sense");
        }
        this.dimension = dimension;
        this.lattice = new HashMap<>(initialCapacity);
    }

    public Lattice(Lattice lattice) {
        this.dimension = lattice.dimension;
        if (dimension < 2 || dimension > 3) {
            throw new IllegalArgumentException("Dimension of less than 2 or more than 3 does not make sense");
        }
        this.lattice = new HashMap<>(lattice.lattice);
        this.energy = lattice.energy;
        this.surfaceSize = lattice.surfaceSize;
        this.plusXBound  = lattice.plusXBound;
        this.minusXBound = lattice.minusXBound;
        this.plusYBound  = lattice.plusYBound;
        this.minusYBound = lattice.minusYBound;
        this.plusZBound  = lattice.plusZBound;
        this.minusZBound = lattice.minusZBound;
    }

    public int getDimension() {
        return dimension;
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
    public boolean containsPoint(Point point) {
        return lattice.containsKey(point);
    }

    /**
     * Returns true if the specified peptide is in the lattice
     * @param peptide
     * @return
     */
    public boolean containsValue(Peptide peptide) {
        return lattice.containsValue(peptide);
    }

    /**
     * Returns the peptide at the specified point, or null
     * @param point
     * @return
     */
    public Peptide get(Point point) {
        return lattice.get(point);
    }

    /**
     * Puts the peptide in the lattice at the specified point.
     * Returns the peptide that previously occupied the point.
     * @param point
     * @param peptide
     * @return
     */
    public void put(Point point, Peptide peptide) {
        if (dimension == 2 && point.z != 0) {
            throw new IllegalArgumentException("2D points cannot have a z-component");
        }
        if (containsPoint(point)) {
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
        for (Direction d : Direction.values(dimension)) {
            Peptide adj = get(point.getAdjacent(d));
            if (adj != null) {
                if (adj.index >= 0) { // this is for use with surface lattice, so that it properly handles surface-perimeter
                    surfaceSize -= 1;
                }
                // if they are not adjoining peptides
                if (adj.index != peptide.index + 1 && adj.index != peptide.index - 1) {
                    energy += peptide.interaction(adj);
                }
                energy -= adj.interaction(Residue.H2O);
            } else {
                surfaceSize += 1;
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
     * Resets the lattice to an empty state
     */
    public void clear() {
        lattice.clear();
        energy = 0;
        surfaceSize = 0;
        plusXBound  = 0;
        minusXBound = 0;
        plusYBound  = 0;
        minusYBound = 0;
        plusZBound  = 0;
        minusZBound = 0;
    }

    /**
     * Returns a set of all occupied points in the lattice
     * @return
     */
    public Set<Point> keySet() {
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
        if (dimension == 2) {
            return 2 * (xRange + yRange);
        } else {
            int zRange = plusZBound - minusZBound + 1;
            return 2 * (xRange * yRange + xRange * zRange + yRange * zRange);
        }
    }

    /**
     * Draws an ASCII visualization of the peptides in the lattice to the console.
     * Also returns the drawing as a list of strings. Currently only draws 2D lattices.
     *
     * TODO: handle drawing 3D lattice better
     */
    public List<String> visualize() {
        List<String> lines = new ArrayList<>();
        for (int k=minusZBound; k<=plusZBound; k++) {
            for (int i=plusYBound; i>=minusYBound; i--) {
                StringBuilder latticeString = new StringBuilder();
                StringBuilder connectionsString = new StringBuilder();
                for (int j=minusXBound; j<=plusXBound; j++) {
                    Peptide p = get(Point.point(j, i, k));
                    if (p != null) {
                        int index = p.index;
                        String residue = p.residue.toString();
                        Point up = Point.point(j, i, k + 1);
                        if (containsPoint(up) && (get(up).index == index + 1 || get(up).index == index - 1)) {
                            residue = residue.replace('(', '{').replace(')', '}');
                        }
                        latticeString.append(residue);
                        Point right = Point.point(j + 1, i, k);
                        if (containsPoint(right) && (get(right).index == index + 1 || get(right).index == index - 1)) {
                            latticeString.append("-");
                        } else {
                            latticeString.append(" ");
                        }
                        Point below = Point.point(j, i - 1, k);
                        if (containsPoint(below) && (get(below).index == index + 1 || get(below).index == index - 1)) {
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
            StringBuilder edge = new StringBuilder();
            for (int j=minusXBound; j<=plusXBound; j++) {
                edge.append("====");
            }
            lines.add(edge.toString());
            System.out.println(edge);
        }
        return lines;
    }
}
