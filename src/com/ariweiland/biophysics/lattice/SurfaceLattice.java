package com.ariweiland.biophysics.lattice;

import com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.peptide.Residue;

import java.util.ArrayList;
import java.util.List;

/**
 * A surface lattice represents a specific type of lattices with a surface
 * at y == 0 composed of a specified residue. Points at y < 0 are invalid
 * points, and will cause exceptions to be thrown.
 * @author Ari Weiland
 */
public class SurfaceLattice extends Lattice {

    private final Residue surface;

    public SurfaceLattice(int dimension, Residue surface) {
        super(dimension);
        if (surface == Residue.H2O) {
            throw new IllegalArgumentException("Surface cannot be water (null)");
        }
        this.surface = surface;
    }

    public SurfaceLattice(int dimension, Residue surface, int initialCapacity) {
        super(dimension, initialCapacity);
        if (surface == Residue.H2O) {
            throw new IllegalArgumentException("Surface cannot be water (null)");
        }
        this.surface = surface;
    }

    public SurfaceLattice(SurfaceLattice lattice) {
        super(lattice);
        this.surface = lattice.surface;
    }

    @Override
    public boolean containsPoint(Point point) {
        if (point.y < 0) {
            throw new IllegalArgumentException("Surface lattices do not have points below y == 0");
        }
        return point.y == 0 || super.containsPoint(point);
    }

    @Override
    public Peptide get(Point point) {
        if (point.y < 0) {
            throw new IllegalArgumentException("Surface lattices do not have points below y == 0");
        } else if (point.y == 0) {
            return new Peptide(-2, surface);
        } else {
            return super.get(point);
        }
    }

    @Override
    public void put(Point point, Peptide peptide) {
        if (point.y < 1) {
            throw new IllegalArgumentException("Cannot put a point on or below the surface (y <= 0)");
        }
        super.put(point, peptide);
    }

    @Override
    public List<String> visualize() {
        List<String> lines = new ArrayList<>();
        for (int k=minusZBound; k<=plusZBound; k++) {
            for (int i=plusYBound; i > 0; i--) {
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
            StringBuilder surface = new StringBuilder();
            StringBuilder base = new StringBuilder();
            for (int j=minusXBound; j<=plusXBound; j++) {
                surface.append(this.surface).append(" ");
                base.append("-+--");
            }
            lines.add(surface.toString());
            lines.add(base.toString());
            System.out.println(surface);
            System.out.println(base);
        }
        return lines;
    }
}
