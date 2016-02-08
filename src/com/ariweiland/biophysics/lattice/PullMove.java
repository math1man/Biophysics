package com.ariweiland.biophysics.lattice;

import com.ariweiland.biophysics.Direction;

/**
 * @author Ari Weiland
 */
public class PullMove {
    public final int index;
    public final Direction direction;

    public PullMove(int index, Direction direction) {
        this.index = index;
        this.direction = direction;
    }

    @Override
    public String toString() {
        return "PullMove{" +
                "index=" + index +
                ", direction=" + direction +
                '}';
    }
}
