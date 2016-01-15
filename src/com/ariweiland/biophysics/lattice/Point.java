package com.ariweiland.biophysics.lattice;

/**
 * @author Ari Weiland
 */
public interface Point {

    public abstract Point getAdjacent(Direction direction);
}
