package com.ariweiland.biophysics;

/**
 * Simple wrapper class for a coordinate in a lattice.
 * Also has a convenience method to get adjacent points.
 * @author Ari Weiland
 */
public class Point {
    public final int x;
    public final int y;

    public Point(int x, int y) {
        this.x = x;
        this.y = y;
    }

    /**
     * Returns a new point that is
     * @param direction
     * @return
     */
    public Point getAdjacent(Direction direction) {
        switch (direction) {
            case EAST:
                return new Point(x + 1, y);
            case NORTH:
                return new Point(x, y + 1);
            case WEST:
                return new Point(x - 1, y);
            case SOUTH:
            default:
                return new Point(x, y - 1);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Point)) return false;

        Point point = (Point) o;

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

    public static enum Direction {
        NORTH, EAST, SOUTH, WEST;

        public Direction getLeft() {
            return values()[(ordinal() + 3) % 4];
        }

        public Direction getReverse() {
            return values()[(ordinal() + 2) % 4];
        }

        public Direction getRight() {
            return values()[(ordinal() + 1) % 4];
        }
    }
}
