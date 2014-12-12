package com.ariweiland.biophysics;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * @author Ari Weiland
 */
public class Lattice {

    private final Map<Point, Peptide> lattice;
    private int plusXBound;
    private int minusXBound;
    private int plusYBound;
    private int minusYBound;
    private double energy = 0;
    private int perimeter = 0;

    public Lattice() {
        lattice = new HashMap<Point, Peptide>();
    }

    public Lattice(Map<Point, Peptide> lattice) {
        this();
        putAll(lattice);
    }

    public Lattice(Lattice lattice) {
        this.lattice = new HashMap<Point, Peptide>(lattice.lattice);
        this.plusXBound = lattice.plusXBound;
        this.minusXBound = lattice.minusXBound;
        this.plusYBound = lattice.plusYBound;
        this.minusYBound = lattice.minusYBound;
        this.energy = lattice.energy;
        this.perimeter = lattice.perimeter;
    }

    public int size() {
        return lattice.size();
    }

    public boolean isEmpty() {
        return lattice.isEmpty();
    }

    public boolean containsPoint(int x, int y) {
        return containsPoint(new Point(x, y));
    }

    public boolean containsPoint(Point point) {
        return lattice.containsKey(point);
    }

    public boolean containsValue(Peptide p) {
        return lattice.containsValue(p);
    }

    public Peptide get(int x, int y) {
        return get(new Point(x, y));
    }

    public Peptide get(Point point) {
        return lattice.get(point);
    }

    public Peptide put(int x, int y, Peptide p) {
        return put(new Point(x, y), p);
    }

    public Peptide put(Point point, Peptide p) {
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
                if (adj.getIndex() != p.getIndex() + 1 && adj.getIndex() != p.getIndex() - 1) {
                    energy += p.interaction(adj);
                }
                energy -= adj.interaction(null);
            } else {
                perimeter += 1;
                energy += p.interaction(null);
            }
        }
        return lattice.put(point, p);
    }

    public void putAll(Map<Point, Peptide> lattice) {
        for (Map.Entry<Point, Peptide> e : lattice.entrySet()) {
            put(e.getKey(), e.getValue());
        }
    }

    public Peptide remove(int x, int y) {
        return remove(new Point(x, y));
    }

    public Peptide remove(Point point) {
        return lattice.remove(point);
    }

    public void clear() {
        lattice.clear();
    }

    public Set<Point> keySet() {
        return lattice.keySet();
    }

    public double getEnergy() {
        return Math.round(energy * 100) / 100.0;
    }

    public int getPerimeter() {
        return perimeter;
    }

    public int getMaxPerim() {
        return 2 * (plusXBound - minusXBound + plusYBound - minusYBound + 2);
    }

    public void visualize() {
        for (int i=minusYBound; i<=plusYBound; i++) {
            String latticeString = "";
            String connectionsString = "";
            for (int j=minusXBound; j<=plusXBound; j++) {
                Peptide p = get(j, i);
                if (p != null) {
                    latticeString += p.getType();
                    int index = p.getIndex();
                    if (containsPoint(j + 1, i) && (get(j + 1, i).getIndex() == index + 1 || get(j + 1, i).getIndex() == index - 1)) {
                        latticeString += "-";
                    } else {
                        latticeString += " ";
                    }
                    if (containsPoint(j, i + 1) && (get(j, i + 1).getIndex() == index + 1 || get(j, i + 1).getIndex() == index - 1)) {
                        connectionsString += " |  ";
                    } else {
                        connectionsString += "    ";
                    }
                } else {
                    latticeString += "    ";
                    connectionsString += "    ";
                }
            }
            System.out.println(latticeString);
            System.out.println(connectionsString);
        }
    }
}
