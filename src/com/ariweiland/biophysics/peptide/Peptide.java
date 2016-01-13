package com.ariweiland.biophysics.peptide;

/**
 * This class wraps a residue with an index, which is useful for a variety of things.
 * It includes wrapper methods for the Residue methods too.
 * @author Ari Weiland
 */
public class Peptide {

    public final int index;
    public final Residue residue;

    public Peptide(int index, Residue residue) {
        this.index = index;
        this.residue = residue;
    }

    public double interaction(Peptide p) {
        return interaction(p == null ? Residue.H2O : p.residue);
    }

    public double interaction(Residue r) {
        return residue.interaction(r);
    }

    public double minInteraction() {
        return residue.minInteraction();
    }

    @Override
    public String toString() {
        return index + ": " + residue;
    }
}
