package com.ariweiland.biophysics;

import java.util.*;

/**
 * @author Ari Weiland
 */
public class Polypeptide {
    
    private final List<PType> polypeptide;
    private final Map<PType, Integer> typeCount = new HashMap<PType, Integer>();

    public Polypeptide() {
        this(new ArrayList<PType>());
    }

    public Polypeptide(List<PType> polypeptide) {
        this.polypeptide = polypeptide;
    }

    public int size() {
        return polypeptide.size();
    }

    public boolean isEmpty() {
        return polypeptide.isEmpty();
    }

    public boolean contains(PType p) {
        return polypeptide.contains(p);
    }

    public Iterator<PType> iterator() {
        return polypeptide.iterator();
    }

    public boolean add(PType type) {
        if (!typeCount.containsKey(type)) {
            typeCount.put(type, 0);
        }
        typeCount.put(type, 1 + typeCount.get(type));
        return polypeptide.add(type);
    }

    public boolean addAll(Collection<? extends PType> c) {
        boolean changed = false;
        for (PType t : c) {
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

    public int indexOf(PType type) {
        return polypeptide.indexOf(type);
    }

    public int lastIndexOf(PType type) {
        return polypeptide.lastIndexOf(type);
    }

    public double getMinEnergy() {
        double minEnergy = 0;
        for (Map.Entry<PType, Integer> e : typeCount.entrySet()) {
            minEnergy += 2 * e.getKey().minInteraction() * e.getValue();
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
