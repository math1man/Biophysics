package com.ariweiland.biophysics;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Ari Weiland
 */
public class PType {

    public static final PType POS = new PType("(+)");
    public static final PType NEG = new PType("(-)");
    public static final PType P = new PType("(P)");
    public static final PType H = new PType("(H)");
    public static final PType NEUT = new PType("( )");
    public static final PType H2O = null;
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

    private PType(String symbol) {
        this.symbol = symbol;
    }

    public double interaction(PType t) {
        return interaction(this, t);
    }

    public double minInteraction() {
        return minInteraction(this);
    }

    @Override
    public String toString() {
        return symbol;
    }

    public static double interaction(PType t1, PType t2) {
        return ENERGY_MAP.get(t1, t2);
    }

    public static double minInteraction(PType t) {
        return ENERGY_MAP.getMin(t);
    }

    private static class EnergyMap {
        private final Map<PType, Map<PType, Double>> map = new HashMap<PType, Map<PType, Double>>();
        
        public void put(PType t1, PType t2, double d) {
            if (!map.containsKey(t1)) {
                map.put(t1, new HashMap<PType, Double>());
            }
            map.get(t1).put(t2, d);
            if (!map.containsKey(t2)) {
                map.put(t2, new HashMap<PType, Double>());
            }
            map.get(t2).put(t1, d);
        }
        
        public double get(PType t1, PType t2) {
            if (map.containsKey(t1)) {
                Map<PType, Double> sub = map.get(t1);
                if (sub.containsKey(t2)) {
                    return sub.get(t2);
                }
            }
            return 0;
        }

        public double getMin(PType t) {
            double min = 0;
            if (map.containsKey(t)) {
                Map<PType, Double> sub = map.get(t);
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
