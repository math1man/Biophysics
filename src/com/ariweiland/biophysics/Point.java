package com.ariweiland.biophysics;

import java.util.Arrays;

/**
 * Simple wrapper class for a coordinate in a lattice.
 * Also has a convenience method to get adjacent points.
 * @author Ari Weiland
 */
public class Point {

    private final int[] coords;

    public Point(int... coords) {
        this.coords = coords;
    }

    public int[] getCoords() {
        return Arrays.copyOf(coords, coords.length);
    }

    /**
     * Returns a new point that is
     * @param direction
     * @return
     */
    public Point getAdjacent(Direction direction) {
        int d = getDimension();
        int[] adj = getCoords();
        if (d > 0 && direction == Direction.EAST) {
            adj[0] += 1;
        } else if (d > 0 && direction == Direction.WEST) {
            adj[0] -= 1;
        } else if (d > 1 && direction == Direction.NORTH) {
            adj[1] += 1;
        } else if (d > 1 && direction == Direction.SOUTH) {
            adj[1] -= 1;
        } else if (d > 2 && direction == Direction.UP) {
            adj[2] += 1;
        } else if (d > 2 && direction == Direction.DOWN) {
            adj[2] -= 1;
        } else {
            throw new IllegalArgumentException(direction + " is out of the dimension of this point");
        }
        return new Point(adj);
    }

    public int getDimension() {
        return getCoords().length;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Point point = (Point) o;

        return Arrays.equals(getCoords(), point.getCoords());

    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(getCoords());
    }

    @Override
    public String toString() {
        return Arrays.toString(getCoords());
    }
}
