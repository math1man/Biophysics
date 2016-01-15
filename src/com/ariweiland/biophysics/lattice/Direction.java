package com.ariweiland.biophysics.lattice;

/**
* @author Ari Weiland
*/
public enum Direction {
    NORTH, EAST, UP, SOUTH, WEST, DOWN;

    public Direction getReverse() {
        return values()[(ordinal() + 3) % 6];
    }

    public static Direction[] values2D() {
        return new Direction[]{ NORTH, EAST, SOUTH, WEST };
    }

    public static Direction[] values3D() {
        return values();
    }
}
