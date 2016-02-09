package com.ariweiland.biophysics.lattice;

import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.peptide.Residue;

import java.util.Stack;

/**
 * @author Ari Weiland
 */
public class BacktrackLattice extends FastLattice {

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
    public void put(Point point, Peptide peptide) {
        energyStack.push(energy);
        super.put(point, peptide);
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
