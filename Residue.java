package com.ariweiland.biophysics;

import java.util.HashMap;
import java.util.Map;

/**
 * This class enumerates the different types of Peptide residues and their energy interactions.
 * @author Ari Weiland
 */
public class Residue {

    public static final Residue POS = new Residue("(+)");
    public static final Residue NEG = new Residue("(-)");
    public static final Residue P = new Residue("(P)");
    public static final Residue H = new Residue("(H)");
    public static final Residue NEUT = new Residue("( )");
    public static final Residue H2O = null;
    private static final EnergyMap ENERGY_MAP = new EnergyMap();

    static {
        // all in ev/kT for T=310K
        // note that favorable (negative) interactions with water disrupt the algorithm
        // ion-ion: +-1.24
        ENERGY_MAP.put(POS, POS,  1.24);
        ENERGY_MAP.put(POS, NEG, -1.24);
        ENERGY_MAP.put(NEG, NEG,  1.24);

        // ion-dipole: ~0.62
        ENERGY_MAP.put(POS, P,   -0.62);
        ENERGY_MAP.put(NEG, P,   -0.62);
//        ENERGY_MAP.put(POS, H2O, -0.62);
//        ENERGY_MAP.put(NEG, H2O, -0.62);

        // dipole-dipole: ~0.11
        ENERGY_MAP.put(P,   P,   -0.11);
//        ENERGY_MAP.put(P,   H2O, -0.11);

        // hydrophobic: +-1.16 (?)
        // (assuming 3 kJ/mol)
        ENERGY_MAP.put(H,   H,   -1.16);
        ENERGY_MAP.put(H,   H2O,  1.16);

        // TODO: add more interactions with H2O?
    }

    private final String symbol;

    private Residue(String symbol) {
        this.symbol = symbol;
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
        
        public void put(Residue t1, Residue t2, double d) {
            if (!map.containsKey(t1)) {
                map.put(t1, new HashMap<Residue, Double>());
            }
            map.get(t1).put(t2, d);
            if (!map.containsKey(t2)) {
                map.put(t2, new HashMap<Residue, Double>());
            }
            map.get(t2).put(t1, d);
        }
        
        public double get(Residue t1, Residue t2) {
            if (map.containsKey(t1)) {
                Map<Residue, Double> sub = map.get(t1);
                if (sub.containsKey(t2)) {
                    return sub.get(t2);
                }
            }
            return 0;
        }

        public double getMin(Residue t) {
            double min = 0;
            if (map.containsKey(t)) {
                Map<Residue, Double> sub = map.get(t);
                min = get(t, t);
                for (double d : sub.values()) {
                    if (d < min) {
                        min = d;
                    }
                }
            }
            return min;
        }
    }
}
