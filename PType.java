package com.ariweiland.biophysics;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Ari Weiland
 */
public class PType {

    public static final PType P = new PType("(P)");
    public static final PType H = new PType("(H)");
    public static final PType H2O = null;
    private static final EnergyMap ENERGY_MAP = new EnergyMap();

    static {
        ENERGY_MAP.put(H, H,  -1.0);
        ENERGY_MAP.put(H, H2O, 1.0);
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
