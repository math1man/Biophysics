package com.ariweiland.biophysics.peptide;

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

    public static final double ION_ION = 1.24;
    public static final double ION_DIPOLE = 0.62;
    public static final double DIPOLE_DIPOLE = 0.62;
    public static final double HYDROPHOBIC = 1; // 1.16

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

    private static double hhInt = 0;
    private static double hpInt = 0;
    private static double ppInt = 0;
    private static double shInt = 0;
    private static double spInt = 0;

    public static void setInteractionScheme(double hh, double hp, double pp) {
        hhInt = hh;
        hpInt = hp;
        ppInt = pp;
    }

    public static void setSurfaceInteractions(double sh, double sp) {
        shInt = sh;
        spInt = sp;
    }

    public static double interaction(Residue r1, Residue r2) {
        if (r1 == P) {
            if (r2 == H) {
                return hpInt;
            } else if (r2 == P) {
                return ppInt;
            } else if (r2 == S) {
                return spInt;
            }
        } else if (r1 == H) {
            if (r2 == H) {
                return hhInt;
            } else if (r2 == P) {
                return hpInt;
            } else if (r2 == S) {
                return shInt;
            }
        }
        return 0;
    }

    private static double min(double... values) {
        double min = values[0];
        for (double v : values) {
            min = Math.min(min, v);
        }
        return min;
    }

    public static double minInteraction(Residue residue) {
        if (residue == P) {
            return min(hpInt, ppInt, spInt, 0);
        } else if (residue == H) {
            return min(hpInt, hhInt, shInt, 0);
        } else {
            return 0;
        }
    }
}
