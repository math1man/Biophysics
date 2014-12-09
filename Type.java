package com.ariweiland.biophysics;

import java.util.HashMap;
import java.util.Map;

/**
* @author Ari Weiland
*/
public enum Type {
    POSITIVE("(+)"), NEGATIVE("(-)"), POLAR("(P)"), NONPOLAR("(N)"), NEUTRAL("( )");

    public static final double CHARGE_INTERACTION_ENERGY = 2;
    public static final double POLAR_INTERACTION_ENERGY = -0;
    public static final double NONPOLAR_INTERACTION_ENERGY = -1;
    public static final double CROSS_POLAR_INTERACTION_ENERGY = 0;
    public static final double H2O_POLAR_INTERACTION_ENERGY = -0;
    public static final double H2O_NONPOLAR_INTERACTION_ENERGY = 1;

    private static final Map<Pairing, Double> ENERGIES = new HashMap<Pairing, Double>();

    static {
        ENERGIES.put(new Pairing(POSITIVE,  POSITIVE),  CHARGE_INTERACTION_ENERGY);
        ENERGIES.put(new Pairing(NEGATIVE,  NEGATIVE),  CHARGE_INTERACTION_ENERGY);
        ENERGIES.put(new Pairing(POSITIVE,  NEGATIVE), -CHARGE_INTERACTION_ENERGY);
        ENERGIES.put(new Pairing(POLAR,     POLAR),     POLAR_INTERACTION_ENERGY);
        ENERGIES.put(new Pairing(POLAR,     NONPOLAR),  CROSS_POLAR_INTERACTION_ENERGY);
        ENERGIES.put(new Pairing(POLAR,     null),      H2O_POLAR_INTERACTION_ENERGY);
        ENERGIES.put(new Pairing(NONPOLAR,  NONPOLAR),  NONPOLAR_INTERACTION_ENERGY);
        ENERGIES.put(new Pairing(NONPOLAR,  null),      H2O_NONPOLAR_INTERACTION_ENERGY);
    }

    private final String symbol;

    Type(String symbol) {
        this.symbol = symbol;
    }

    @Override
    public String toString() {
        return symbol;
    }

    public double interaction(Type p) {
        return interaction(this, p);
    }

    public double minInteraction() {
        double minEnergy = interaction(null);
        for (Type t : values()) {
            double energy = interaction(t);
            if (energy < minEnergy) {
                minEnergy = energy;
            }
        }
        return minEnergy;
    }

    public static double interaction(Type p1, Type p2) {
        Pairing pairing = new Pairing(p1, p2);
        return ENERGIES.containsKey(pairing) ? ENERGIES.get(pairing) : 0;
    }

    private static class Pairing {
        final Type t1;
        final Type t2;

        private Pairing(Type t1, Type t2) {
            this.t1 = t1;
            this.t2 = t2;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof Pairing)) return false;

            Pairing pairing = (Pairing) o;

            return t1 == pairing.t1 && t2 == pairing.t2 || t1 == pairing.t2 && t2 == pairing.t1;

        }

        @Override
        public int hashCode() {
            // needs to be symmetrical
            return (t1 != null ? t1.hashCode() : 0) + (t2 != null ? t2.hashCode() : 0);
        }
    }
}
