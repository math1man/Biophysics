package com.ariweiland.biophysics;

/**
 * Simple wrapper class for a coordinate in a lattice.
 * Also has a convenience method to get adjacent points.
 * @author Ari Weiland
 */
public class Point {

    public final int x;
    public final int y;
    public final int z;

    public Point(int x, int y) {
        this(x, y, 0);
    }

    public Point(int x, int y, int z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    /**
     * Returns a new point that is
     * @param direction
     * @return
     */
    public Point getAdjacent(Direction direction) {
        switch (direction) {
            case EAST:
                return new Point(x + 1, y, z);
            case WEST:
                return new Point(x - 1, y, z);
            case NORTH:
                return new Point(x, y + 1, z);
            case SOUTH:
                return new Point(x, y - 1, z);
            case UP:
                return new Point(x, y, z + 1);
            case DOWN:
                return new Point(x, y, z - 1);
            default: // should never get called
                return this;
        }
    }
    
    public boolean isAdjacentTo(Point p) {
        int shifts = 0;
        int dif = Math.abs(x - p.x);
        if (dif == 1) {
            shifts++;
        } else if (dif > 1) {
            return false;
        }
        dif = Math.abs(y - p.y);
        if (dif == 1) {
            shifts++;
        } else if (dif > 1) {
            return false;
        }
        dif = Math.abs(z - p.z);
        if (dif == 1) {
            shifts++;
        } else if (dif > 1) {
            return false;
        }
        return shifts == 1;
    }

    public Direction getDirectionTo(Point p) {
        if (!isAdjacentTo(p)) {
            throw new IllegalArgumentException("Points are not adjacent");
        }
        if (p.x > x) {
            return Direction.EAST;
        } else if (p.x < x) {
            return Direction.WEST;
        } else if (p.y > y) {
            return Direction.NORTH;
        } else if (p.y < y) {
            return Direction.SOUTH;
        } else if (p.z > z) {
            return Direction.UP;
        } else {
            return Direction.DOWN;
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Point point = (Point) o;

        return x == point.x && y == point.y && z == point.z;

    }

    @Override
    public int hashCode() {
        int hashCode = (x >>> 20) ^ (x >>> 12) ^ x ^ (y >>> 22) ^ (y >>> 10) ^ y ^ (z >>> 18) ^ (z >>> 14) ^ z;
        return hashCode ^ (hashCode >>> 7) ^ (hashCode >>> 4);
    }

    @Override
    public String toString() {
        return "(" + x + ", " + y + ", " + z + ")";
    }
}
