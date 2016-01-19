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
    protected final int plusBounds[];
    protected final int minusBounds[];

    public Lattice(int dimension) {
        if (dimension < 2 || dimension > 3) {
            throw new IllegalArgumentException("Dimension of less than 2 or more than 3 does not make sense");
        }
        this.dimension = dimension;
        this.lattice = new HashMap<>();
        plusBounds = new int[dimension];
        minusBounds = new int[dimension];
    }

    public Lattice(Lattice lattice) {
        this.dimension = lattice.dimension;
        if (dimension < 2 || dimension > 3) {
            throw new IllegalArgumentException("Dimension of less than 2 or more than 3 does not make sense");
        }
        this.lattice = new HashMap<>(lattice.lattice);
        this.energy = lattice.energy;
        this.plusBounds = Arrays.copyOf(lattice.plusBounds, dimension);
        this.minusBounds = Arrays.copyOf(lattice.minusBounds, dimension);
        this.surfaceSize = lattice.surfaceSize;
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
        if (point.getDimension() != dimension) {
            throw new IllegalArgumentException("Incorrect dimension");
        }
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
        if (point.getDimension() != dimension) {
            throw new IllegalArgumentException("Incorrect number of points");
        }
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
        if (containsPoint(point)) {
            throw new IllegalArgumentException("That point is already occupied");
        }
        if (isEmpty()) {
            System.arraycopy(point.getCoords(), 0, plusBounds, 0, dimension);
            System.arraycopy(point.getCoords(), 0, minusBounds, 0, dimension);
        } else {
            for (int i=0; i<dimension; i++) {
                if (point.getCoords()[i] > plusBounds[i]) {
                    plusBounds[i] = point.getCoords()[i];
                } else if (point.getCoords()[i] < minusBounds[i]) {
                    minusBounds[i] = point.getCoords()[i];
                }
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
        Arrays.fill(plusBounds, 0);
        Arrays.fill(minusBounds, 0);
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
        int bound = 0;
        for (int i=0; i<dimension; i++) {
            bound += plusBounds[i] - minusBounds[i] + 1;
        }
        return 2 * bound;
    }

    /**
     * Draws an ASCII visualization of the peptides in the lattice to the console.
     * Also returns the drawing as a list of strings. Currently only draws 2D lattices.
     *
     * TODO: handle drawing 3D lattice???
     */
    public List<String> visualize() {
        List<String> lines = new ArrayList<>();
        if (dimension == 2) {
            for (int i=plusBounds[1]; i>=minusBounds[1]; i--) {
                StringBuilder latticeString = new StringBuilder();
                StringBuilder connectionsString = new StringBuilder();
                for (int j=minusBounds[0]; j<=plusBounds[0]; j++) {
                    Peptide p = get(new Point(j, i));
                    if (p != null) {
                        latticeString.append(p.residue);
                        int index = p.index;
                        Point right = new Point(j + 1, i);
                        if (containsPoint(right) && (get(right).index == index + 1 || get(right).index == index - 1)) {
                            latticeString.append("-");
                        } else {
                            latticeString.append(" ");
                        }
                        Point below = new Point(j, i - 1);
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
        }
        return lines;
    }
}
