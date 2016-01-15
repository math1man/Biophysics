package com.ariweiland.biophysics;

/**
* @author Ari Weiland
*/
public enum Direction {
    NORTH, EAST, SOUTH, WEST;

    public Direction getReverse() {
        return values()[(ordinal() + 2) % 4];
    }

}
