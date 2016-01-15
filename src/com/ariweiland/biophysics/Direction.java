package com.ariweiland.biophysics;

/**
* @author Ari Weiland
*/
public enum Direction {
    EAST, WEST, NORTH, SOUTH, UP, DOWN;

    public Direction getReverse() {
        return values()[(ordinal() + (ordinal() % 2 == 0 ? 1 : -1))];
    }

    public static Direction[] values(int dimensions) {
        Direction[] directions = new Direction[2 * dimensions];
        System.arraycopy(values(), 0, directions, 0, 2 * dimensions);
        return directions;
    }

}
