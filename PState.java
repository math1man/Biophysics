package com.ariweiland.biophysics;

/**
* @author Ari Weiland
*/
class PState implements Comparable<PState> {
    final Lattice lattice;
    final Point point;
    final int index;
    final double energyBound;
    final int perimeter;

    PState(Lattice lattice, int x, int y, int index, double energyBound) {
        this(lattice, new Point(x, y), index, energyBound);
    }

    PState(Lattice lattice, Point point, int index, double energyBound) {
        this.lattice = lattice;
        this.point = point;
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
