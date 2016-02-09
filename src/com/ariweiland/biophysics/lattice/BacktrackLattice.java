package com.ariweiland.biophysics.lattice;

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
        this(dimension, null);
    }

    public BacktrackLattice(int dimension, Residue surface) {
        super(dimension, surface);
    }

    public BacktrackLattice(int dimension, int initialCapacity) {
        this(dimension, initialCapacity, null);
    }

    public BacktrackLattice(int dimension, int initialCapacity, Residue surface) {
        super(dimension, initialCapacity, surface);
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
