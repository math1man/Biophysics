package com.ariweiland.biophysics;

/**
 * @author Ari Weiland
 */
public class Peptide {

    private final int index;
    private final PType type;

    public Peptide(int index, PType type) {
        this.index = index;
        this.type = type;
    }

    public int getIndex() {
        return index;
    }

    public PType getType() {
        return type;
    }

    public double interaction(Peptide p) {
        return type.interaction(p == null ? PType.H2O : p.getType());
    }

    public double minInteraction() {
        return type.minInteraction();
    }

    @Override
    public String toString() {
        return index + ": " + type;
    }
}
