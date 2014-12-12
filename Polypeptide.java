package com.ariweiland.biophysics;

import java.util.*;

/**
 * @author Ari Weiland
 */
public class Polypeptide {

    public static final Polypeptide GLUCAGON = new Polypeptide("+PP PHPP-HP+HH-P++HP-HHPHHHPP");
    
    private final List<PType> polypeptide;
    private final Map<PType, Integer> typeCount = new HashMap<PType, Integer>();

    public Polypeptide() {
        this(new ArrayList<PType>());
    }

    public Polypeptide(List<PType> polypeptide) {
        this.polypeptide = polypeptide;
    }

    public Polypeptide(String peptideString) {
        this();
        peptideString = peptideString.toUpperCase();
        for (int i=0; i<peptideString.length(); i++) {
            switch (peptideString.charAt(i)) {
                case '+':
                    add(PType.POS);
                    break;
                case '-':
                    add(PType.NEG);
                    break;
                case 'P':
                    add(PType.P);
                    break;
                case 'H':
                    add(PType.H);
                    break;
                case ' ':
                    add(PType.NEUT);
                    break;
                default:
                    break;
            }
        }
    }

    public int size() {
        return polypeptide.size();
    }

    public boolean isEmpty() {
        return polypeptide.isEmpty();
    }

    public boolean add(PType type) {
        if (!typeCount.containsKey(type)) {
            typeCount.put(type, 0);
        }
        typeCount.put(type, 1 + typeCount.get(type));
        return polypeptide.add(type);
    }

    public void clear() {
        polypeptide.clear();
    }

    public Peptide get(int index) {
        return new Peptide(index, polypeptide.get(index));
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
