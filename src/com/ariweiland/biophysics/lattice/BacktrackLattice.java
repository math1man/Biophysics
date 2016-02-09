package com.ariweiland.biophysics.lattice;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.peptide.Residue;

import java.util.Stack;

/**
 * @author Ari Weiland
 */
public class BacktrackLattice extends Lattice {

    private final Stack<Double> energyStack = new Stack<>();
    private final Stack<Point> pointStack = new Stack<>();

    public BacktrackLattice(int dimension) {
        super(dimension);
    }

    public BacktrackLattice(int dimension, int initialCapacity) {
        super(dimension, initialCapacity);
    }

    public BacktrackLattice(int dimension, Residue surface) {
        super(dimension, surface);
    }

    public BacktrackLattice(int dimension, Residue surface, int initialCapacity) {
        super(dimension, surface, initialCapacity);
    }

    @Override
    public Peptide get(Point point) {
        if (hasSurface() && point.y == 0) {
            return new Peptide(-2, getSurface());
        } else {
            return lattice.get(point);
        }
    }

    @Override
    public void put(Point point, Peptide peptide) {
        energyStack.push(energy);
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
        pointStack.push(point);
    }

    public void removeLast() {
        energy = energyStack.pop();
        lattice.remove(pointStack.pop());
    }

    public Point getLastPoint() {
        return pointStack.peek();
    }
}
