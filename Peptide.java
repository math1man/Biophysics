package com.ariweiland.biophysics;

/**
 * @author Ari Weiland
 */
public class Peptide {
    private final int index;
    private final Type type;

    public Peptide(int index, Type type) {
        this.index = index;
        this.type = type;
    }

    public int getIndex() {
        return index;
    }

    public Type getType() {
        return type;
    }

    public double interaction(Peptide p) {
        return type.interaction(p == null ? null : p.getType());
    }

    public double minInteraction() {
        return type.minInteraction();
    }

    @Override
    public String toString() {
        return index + ": " + type;
    }
}
