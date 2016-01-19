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
    public final int dimension;

    public Point(Point p) {
        this(p.x, p.y, p.z, p.dimension);
    }

    public Point(int x, int y) {
        this(x, y, 0, 2);
    }

    public Point(int x, int y, int z) {
        this(x, y, z, 3);
    }

    private Point(int x, int y, int z, int dimension) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.dimension = dimension;
    }

    /**
     * Returns a new point that is
     * @param direction
     * @return
     */
    public Point getAdjacent(Direction direction) {
        if (direction == Direction.EAST) {
            return new Point(x + 1, y, z, dimension);
        } else if (direction == Direction.WEST) {
            return new Point(x - 1, y, z, dimension);
        } else if (direction == Direction.NORTH) {
            return new Point(x, y + 1, z, dimension);
        } else if (direction == Direction.SOUTH) {
            return new Point(x, y - 1, z, dimension);
        } else if (dimension > 2 && direction == Direction.UP) {
            return new Point(x, y, z + 1, dimension);
        } else if (dimension > 2 && direction == Direction.DOWN) {
            return new Point(x, y, z - 1, dimension);
        } else {
            throw new IllegalArgumentException(direction + " is out of the dimension of this point");
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Point point = (Point) o;

        if (dimension != point.dimension) return false;
        if (x != point.x) return false;
        if (y != point.y) return false;
        if (z != point.z) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = x;
        result = 31 * result + y;
        result = 31 * result + z;
        result = 31 * result + dimension;
        return result;
    }

    @Override
    public String toString() {
        if (dimension == 2) {
            return "(" + x + ", " + y + ")";
        } else {
            return "(" + x + ", " + y + ", " + z + ")";
        }
    }
}
