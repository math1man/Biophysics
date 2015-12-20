package com.ariweiland.biophysics;

/**
 * This class wraps a residue with an index, which is useful for a variety of things.
 * It includes wrapper methods for the Residue methods too.
 * @author Ari Weiland
 */
public class Peptide {

    public final int index;
    public final Residue type;

    public Peptide(int index, Residue type) {
        this.index = index;
        this.type = type;
    }

    public double interaction(Peptide p) {
        return type.interaction(p == null ? Residue.H2O : p.type);
    }

    public double minInteraction() {
        return type.minInteraction();
    }

    @Override
    public String toString() {
        return index + ": " + type;
    }
}
