package com.ariweiland.biophysics;

import java.util.*;

/**
 * This class represents a folded polypeptide, laid out in 2D rectangular grid.
 * It also keeps track of various relevant properties, such as the bounding box
 * and perimeter of the polypeptide in the lattice and the lattice energy.
 * @author Ari Weiland
 */
public class Lattice {

    private final Map<Point, Peptide> lattice;
    private int plusXBound = 0;
    private int minusXBound = 0;
    private int plusYBound = 0;
    private int minusYBound = 0;
    private int perimeter = 0;
    private double energy = 0;

    public Lattice() {
        this(new HashMap<Point, Peptide>());
    }

    public Lattice(Map<Point, Peptide> lattice) {
        this.lattice = new HashMap<>();
        putAll(lattice);
    }

    public Lattice(Lattice lattice) {
        this.lattice = new HashMap<>(lattice.lattice);
        this.plusXBound = lattice.plusXBound;
        this.minusXBound = lattice.minusXBound;
        this.plusYBound = lattice.plusYBound;
        this.minusYBound = lattice.minusYBound;
        this.perimeter = lattice.perimeter;
        this.energy = lattice.energy;
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
     * @param x
     * @param y
     * @return
     */
    public boolean containsPoint(int x, int y) {
        return containsPoint(new Point(x, y));
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
     * @param x
     * @param y
     * @return
     */
    public Peptide get(int x, int y) {
        return get(new Point(x, y));
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
     * @param x
     * @param y
     * @param peptide
     * @return
     */
    public Peptide put(int x, int y, Peptide peptide) {
        return put(new Point(x, y), peptide);
    }

    /**
     * Puts the peptide in the lattice at the specified point.
     * Returns the peptide that previously occupied the point.
     * @param point
     * @param peptide
     * @return
     */
    public Peptide put(Point point, Peptide peptide) {
        if (isEmpty()) {
            plusXBound = point.x;
            minusXBound = point.x;
            plusYBound = point.y;
            minusYBound = point.y;
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
        }
        for (Point.Direction d : Point.Direction.values()) {
            Peptide adj = get(point.getAdjacent(d));
            if (adj != null) {
                perimeter -= 1;
                if (adj.index != peptide.index + 1 && adj.index != peptide.index - 1) {
                    energy += peptide.interaction(adj);
                }
                energy -= adj.interaction(null);
            } else {
                perimeter += 1;
                energy += peptide.interaction(null);
            }
        }
        return lattice.put(point, peptide);
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
        plusXBound = 0;
        minusXBound = 0;
        plusYBound = 0;
        minusYBound = 0;
        perimeter = 0;
        energy = 0;
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
     * Returns the perimeter of the peptides in the lattice
     * @return
     */
    public int getPerimeter() {
        return perimeter;
    }

    /**
     * Returns the perimeter of the smallest bounding
     * box that contains the peptides in the lattice
     * @return
     */
    public int boundingPerimeter() {
        return 2 * (plusXBound - minusXBound + plusYBound - minusYBound + 2);
    }

    /**
     * Draws an ASCII visualization of the peptides in the lattice to the console.
     * Also returns the drawing as a list of strings.
     */
    public List<String> visualize() {
        List<String> lines = new ArrayList<>();
        for (int i=minusYBound; i<=plusYBound; i++) {
            String latticeString = "";
            String connectionsString = "";
            for (int j=minusXBound; j<=plusXBound; j++) {
                Peptide p = get(j, i);
                if (p != null) {
                    latticeString += p.type;
                    int index = p.index;
                    if (containsPoint(j + 1, i) && (get(j + 1, i).index == index + 1 || get(j + 1, i).index == index - 1)) {
                        latticeString += "-";
                    } else {
                        latticeString += " ";
                    }
                    if (containsPoint(j, i + 1) && (get(j, i + 1).index == index + 1 || get(j, i + 1).index == index - 1)) {
                        connectionsString += " |  ";
                    } else {
                        connectionsString += "    ";
                    }
                } else {
                    latticeString += "    ";
                    connectionsString += "    ";
                }
            }
            lines.add(latticeString);
            lines.add(connectionsString);
            System.out.println(latticeString);
            System.out.println(connectionsString);
        }
        return lines;
    }
}
