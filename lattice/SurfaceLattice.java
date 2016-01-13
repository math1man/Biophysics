package com.ariweiland.biophysics.lattice;

import com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.peptide.Residue;

import java.util.List;

/**
 * A surface lattice represents a specific type of lattices with a surface
 * at y == 0 composed of a specified residue. Points at y < 0 are invalid
 * points, and will cause exceptions to be thrown.
 * @author Ari Weiland
 */
public class SurfaceLattice extends Lattice {

    private final Residue surface;

    public SurfaceLattice(Residue surface) {
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
        List<String> lines = super.visualize();
        for (int i=1; i<minusYBound; i++) {
            lines.add("");
            lines.add("");
            System.out.println();
            System.out.println();
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
        return lines;
    }
}
