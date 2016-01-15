package com.ariweiland.biophysics.lattice;

/**
 * Simple wrapper class for a coordinate in a lattice.
 * Also has a convenience method to get adjacent points.
 * @author Ari Weiland
 */
public class Point2D {
    public final int x;
    public final int y;

    public Point2D(int x, int y) {
        this.x = x;
        this.y = y;
    }

    /**
     * Returns a new point that is
     * @param direction
     * @return
     */
    public Point2D getAdjacent(Direction direction) {
        switch (direction) {
            case EAST:
                return new Point2D(x + 1, y);
            case NORTH:
                return new Point2D(x, y - 1);
            case WEST:
                return new Point2D(x - 1, y);
            case SOUTH:
            default:
                return new Point2D(x, y + 1);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Point2D)) return false;

        Point2D point = (Point2D) o;

        return x == point.x && y == point.y;

    }

    @Override
    public int hashCode() {
        int result = x;
        result = 31 * result + y;
        return result;
    }

    @Override
    public String toString() {
        return "(" + x + ", " + y + ")";
    }

}
