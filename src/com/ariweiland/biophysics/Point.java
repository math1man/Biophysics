package com.ariweiland.biophysics;

import java.util.HashMap;
import java.util.Map;

/**
 * Simple wrapper class for a coordinate in a lattice.
 * Also has a convenience method to get adjacent points.
 * @author Ari Weiland
 */
public class Point {

    public final int x;
    public final int y;
    public final int z;

    private static Map<Integer, Map<Integer, Map<Integer, Point>>> pointMap = new HashMap<>();

    public static Point point(int x, int y) {
        return point(x, y, 0);
    }

    public static Point point(int x, int y, int z) {
        if (!pointMap.containsKey(x)) {
            pointMap.put(x, new HashMap<Integer, Map<Integer, Point>>());
        }
        Map<Integer, Map<Integer, Point>> xMap = pointMap.get(x);
        if (!xMap.containsKey(y)) {
            xMap.put(y, new HashMap<Integer, Point>());
        }
        Map<Integer, Point> yMap = xMap.get(y);
        if (!yMap.containsKey(z)) {
            yMap.put(z, new Point(x, y, z));
        }
        return yMap.get(z);
    }

    private Point(int x, int y, int z) {
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
                return Point.point(x + 1, y, z);
            case WEST:
                return Point.point(x - 1, y, z);
            case NORTH:
                return Point.point(x, y + 1, z);
            case SOUTH:
                return Point.point(x, y - 1, z);
            case UP:
                return Point.point(x, y, z + 1);
            case DOWN:
                return Point.point(x, y, z - 1);
            default: // should never get called
                return this;
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
        int result = x;
        result = 31 * result + y;
        result = 31 * result + z;
        return result;
    }

    @Override
    public String toString() {
        return "(" + x + ", " + y + ", " + z + ")";
    }
}
