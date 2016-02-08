package com.ariweiland.biophysics.peptide;

import java.util.*;

/**
 * This class represents a polypeptide sequence
 * @author Ari Weiland
 */
public class Polypeptide {

    public static final Polypeptide GLUCAGON = new Polypeptide("+PP_PHPPmHP+HHmP++HPmHHPHHHPP");
    
    private final List<Residue> polypeptide;
    private final Map<Residue, Integer> typeCount = new HashMap<>();

    public Polypeptide() {
        this(new ArrayList<Residue>());
    }

    public Polypeptide(List<Residue> polypeptide) {
        this.polypeptide = polypeptide;
    }

    public Polypeptide(String peptideString) {
        this();
        peptideString = peptideString.toUpperCase().replaceAll("[^+mPH_]", "");
        for (int i=0; i<peptideString.length(); i++) {
            add(Residue.get(peptideString.charAt(i)));
        }
    }

    /**
     * Returns the number of peptides that compose this polypeptide
     * @return
     */
    public int size() {
        return polypeptide.size();
    }

    /**
     * Returns true if this polypeptide is empty
     * @return
     */
    public boolean isEmpty() {
        return polypeptide.isEmpty();
    }

    /**
     * Adds a residue to this polypeptide
     * @param type
     */
    public void add(Residue type) {
        if (!typeCount.containsKey(type)) {
            typeCount.put(type, 0);
        }
        typeCount.put(type, 1 + typeCount.get(type));
        polypeptide.add(type);
    }

    /**
     * Removes all peptides from this polypeptide
     */
    public void clear() {
        polypeptide.clear();
    }

    /**
     * Returns the peptide at the given index
     * @param index
     * @return
     */
    public Peptide get(int index) {
        return new Peptide(index, polypeptide.get(index));
    }

    /**
     * Returns the absolute minimum energy of this polypeptide
     * @return
     * @param dimension
     */
    public double getMinEnergy(int dimension) {
        double minEnergy = 0;
        for (Map.Entry<Residue, Integer> e : typeCount.entrySet()) {
            minEnergy += (dimension - 1) * 2 * e.getKey().minInteraction() * e.getValue();
        }
        return minEnergy;
    }

    @Override
    public String toString() {
        String output = "" + polypeptide.get(0);
        for (int i=1; i<size(); i++) {
            output += "-";
            output += polypeptide.get(i);
        }
        return output;
    }
}
