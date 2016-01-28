package com.ariweiland.biophysics;

import java.util.HashMap;
import java.util.Map;

/**
* @author Ari Weiland
*/
public enum Direction {
    EAST, WEST, NORTH, SOUTH, UP, DOWN;

    public Direction getReverse() {
        return values()[(ordinal() + (ordinal() % 2 == 0 ? 1 : -1))];
    }

    public Direction[] getNormals(int dimension) {
        if (dimension == 2) {
            return normals2D.get(this);
        } else {
            return normals3D.get(this);
        }
    }

    public static Direction[] values(int dimensions) {
        Direction[] directions = new Direction[2 * dimensions];
        System.arraycopy(values(), 0, directions, 0, 2 * dimensions);
        return directions;
    }

    private static Map<Direction, Direction[]> normals2D = new HashMap<>();
    private static Map<Direction, Direction[]> normals3D = new HashMap<>();

    static {
        Direction[] ew2D = new Direction[]{NORTH, SOUTH};
        Direction[] ns2D = new Direction[]{EAST, WEST};
        normals2D.put(EAST, ew2D);
        normals2D.put(WEST, ew2D);
        normals2D.put(NORTH, ns2D);
        normals2D.put(SOUTH, ns2D);

        Direction[] ew3D = new Direction[]{NORTH, SOUTH, UP, DOWN};
        Direction[] ns3D = new Direction[]{EAST, WEST, UP, DOWN};
        Direction[] ud3D = new Direction[]{EAST, WEST, NORTH, SOUTH};
        normals3D.put(EAST, ew3D);
        normals3D.put(WEST, ew3D);
        normals3D.put(NORTH, ns3D);
        normals3D.put(SOUTH, ns3D);
        normals3D.put(UP, ud3D);
        normals3D.put(DOWN, ud3D);
    }

}
