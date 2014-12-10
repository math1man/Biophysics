package com.ariweiland.biophysics;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Ari Weiland
 */
public class Type {

    public static final Type P = new Type("(P)");
    public static final Type H = new Type("(H)");
    public static final Type H2O = null;
    private static final EnergyMap ENERGY_MAP = new EnergyMap();

    static {
        ENERGY_MAP.put(H, H,  -1.0);
        ENERGY_MAP.put(H, H2O, 1.0);
    }

    private final String symbol;

    private Type(String symbol) {
        this.symbol = symbol;
    }

    public double interaction(Type t) {
        return interaction(this, t);
    }

    public double minInteraction() {
        return minInteraction(this);
    }

    @Override
    public String toString() {
        return symbol;
    }

    public static double interaction(Type t1, Type t2) {
        return ENERGY_MAP.get(t1, t2);
    }

    public static double minInteraction(Type t) {
        return ENERGY_MAP.getMin(t);
    }

    private static class EnergyMap {
        private final Map<Type, Map<Type, Double>> map = new HashMap<Type, Map<Type, Double>>();
        
        public void put(Type t1, Type t2, double d) {
            if (!map.containsKey(t1)) {
                map.put(t1, new HashMap<Type, Double>());
            }
            map.get(t1).put(t2, d);
            if (!map.containsKey(t2)) {
                map.put(t2, new HashMap<Type, Double>());
            }
            map.get(t2).put(t1, d);
        }
        
        public double get(Type t1, Type t2) {
            if (map.containsKey(t1)) {
                Map<Type, Double> sub = map.get(t1);
                if (sub.containsKey(t2)) {
                    return sub.get(t2);
                }
            }
            return 0;
        }

        public double getMin(Type t) {
            double min = 0;
            if (map.containsKey(t)) {
                Map<Type, Double> sub = map.get(t);
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
