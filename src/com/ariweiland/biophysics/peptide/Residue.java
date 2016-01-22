package com.ariweiland.biophysics.peptide;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * This class enumerates the different types of Peptide residues and their energy interactions.
 * @author Ari Weiland
 */
public class Residue {

    public static final Residue POS = new Residue("(+)");
    public static final Residue NEG = new Residue("(m)");
    public static final Residue P = new Residue("(P)");
    public static final Residue H = new Residue("(H)");
    public static final Residue NEUT = new Residue("(_)");
    public static final Residue S = new Residue("(S)"); // arbitrary "Surface" residue with custom interactions
    public static final Residue H2O = null;
    private static final EnergyMap ENERGY_MAP = new EnergyMap();

    static {
        // all in ev/kT for T=310K
        // ion-ion: +-1.24
        ENERGY_MAP.put(POS, POS,  1.24);
        ENERGY_MAP.put(POS, NEG, -1.24);
        ENERGY_MAP.put(NEG, NEG,  1.24);

        // ion-dipole: ~0.62
        ENERGY_MAP.put(POS, P,   -0.62);
        ENERGY_MAP.put(NEG, P,   -0.62);
        ENERGY_MAP.put(POS, H2O, -0.62);
        ENERGY_MAP.put(NEG, H2O, -0.62);

        // dipole-dipole: ~0.11
//        ENERGY_MAP.put(P,   P,   -0.11);
//        ENERGY_MAP.put(P,   H2O, -0.11);

        // hydrophobic: +-1.16 (?)
        // (assuming 3 kJ/mol)
        ENERGY_MAP.put(H,   H,   -1);
//        ENERGY_MAP.put(H,   H2O,  1.16);
//        ENERGY_MAP.put(H,   P,    1.16);

        ENERGY_MAP.put(S,   H,   -1);
//        ENERGY_MAP.put(S,   P,   -1);

        // TODO: adjust interactions?
    }

    private final String symbol;

    private Residue(String symbol) {
        this.symbol = symbol;
    }

    /**
     * Returns the residue corresponding to the input character.
     * Any characters besides +, m, P, or H default to the neutral residue.
     * @param c
     * @return
     */
    public static Residue get(char c) {
        c = Character.toUpperCase(c);
        switch (c) {
            case '+':
                return Residue.POS;
            case 'M':
                return Residue.NEG;
            case 'P':
                return Residue.P;
            case 'H':
                return Residue.H;
            case 'S':
                return Residue.S;
            default:
                return Residue.NEUT;
        }
    }

    /**
     * Returns the interaction energy between this and another residue, in eV/kT for T=310K
     * @param residue
     * @return
     */
    public double interaction(Residue residue) {
        return interaction(this, residue);
    }

    /**
     * Returns the lowest interaction energy between this and any other residue, in eV/kT for T=310K
     * @return
     */
    public double minInteraction() {
        return minInteraction(this);
    }

    @Override
    public String toString() {
        return symbol;
    }

    public static double interaction(Residue r1, Residue r2) {
        return ENERGY_MAP.get(r1, r2);
    }

    public static double minInteraction(Residue residue) {
        return ENERGY_MAP.getMin(residue);
    }

    private static class EnergyMap {
        private final Map<Residue, Map<Residue, Double>> map = new HashMap<>();
        
        public void put(Residue r1, Residue r2, double d) {
            if (!map.containsKey(r1)) {
                map.put(r1, new HashMap<Residue, Double>());
            }
            map.get(r1).put(r2, d);
            if (r1 != r2) {
                if (!map.containsKey(r2)) {
                    map.put(r2, new HashMap<Residue, Double>());
                }
                map.get(r2).put(r1, d);
            }
        }
        
        public double get(Residue r1, Residue r2) {
            if (map.containsKey(r1)) {
                Map<Residue, Double> sub = map.get(r1);
                if (sub.containsKey(r2)) {
                    return sub.get(r2);
                }
            }
            return 0;
        }

        public double getMin(Residue r) {
            if (map.containsKey(r)) {
                Map<Residue, Double> sub = map.get(r);
                double min = Collections.min(sub.values());
                if (sub.size() < 6) {
                    return Math.min(min, 0);
                } else {
                    return min;
                }
            } else {
                return 0;
            }
        }
    }
}
