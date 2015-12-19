package com.ariweiland.biophysics;

/**
* @author Ari Weiland
*/
class PState implements Comparable<PState> {
    final Lattice lattice;
    final Point lastPoint;
    final int index;
    final double energyBound;
    final int perimeter;

    PState(Lattice lattice, int lastX, int lastY, int index, double energyBound) {
        this(lattice, new Point(lastX, lastY), index, energyBound);
    }

    PState(Lattice lattice, Point lastPoint, int index, double energyBound) {
        this.lattice = lattice;
        this.lastPoint = lastPoint;
        this.index = index;
        this.energyBound = energyBound;
        this.perimeter = lattice.getPerimeter();
    }

    @Override
    public int compareTo(PState o) {
        int compare = Double.compare(energyBound, o.energyBound);
        if (compare == 0) {
            compare = Integer.compare(perimeter, o.perimeter);
        }
        return compare;
    }

    @Override
    public String toString() {
        return energyBound + "/" + perimeter;
    }
}
