package com.ariweiland.biophysics;

import java.util.*;

/**
 * @author Ari Weiland
 */
public class Polypeptide {
    
    private final List<Type> polypeptide;
    private final Map<Type, Integer> typeCount = new HashMap<Type, Integer>();

    {
        for (Type t : Type.values()) {
            typeCount.put(t, 0);
        }
    }

    public Polypeptide() {
        this(new ArrayList<Type>());
    }

    public Polypeptide(List<Type> polypeptide) {
        this.polypeptide = polypeptide;
    }

    public int size() {
        return polypeptide.size();
    }

    public boolean isEmpty() {
        return polypeptide.isEmpty();
    }

    public boolean contains(Type p) {
        return polypeptide.contains(p);
    }

    public Iterator<Type> iterator() {
        return polypeptide.iterator();
    }

    public boolean add(Type type) {
        typeCount.put(type, 1 + typeCount.get(type));
        return polypeptide.add(type);
    }

    public boolean addAll(Collection<? extends Type> c) {
        boolean changed = false;
        for (Type t : c) {
            changed = changed || polypeptide.add(t);
        }
        return changed;
    }

    public void clear() {
        polypeptide.clear();
    }

    public Peptide get(int index) {
        return new Peptide(index, polypeptide.get(index));
    }

    public int indexOf(Type type) {
        return polypeptide.indexOf(type);
    }

    public int lastIndexOf(Type type) {
        return polypeptide.lastIndexOf(type);
    }

    public double getMinEnergy() {
        double minEnergy = 0;
        minEnergy += 2 * typeCount.get(Type.POSITIVE) * Type.POSITIVE.minInteraction();
        minEnergy += 2 * typeCount.get(Type.NEGATIVE) * Type.NEGATIVE.minInteraction();
        minEnergy += 2 * typeCount.get(Type.POLAR)    * Type.POLAR.minInteraction();
        minEnergy += 2 * typeCount.get(Type.NONPOLAR) * Type.NONPOLAR.minInteraction();
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
