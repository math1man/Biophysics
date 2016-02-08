package com.ariweiland.biophysics.modeler;

import com.ariweiland.biophysics.RandomUtils;
import com.ariweiland.biophysics.lattice.Folding;
import com.ariweiland.biophysics.lattice.Lattice;
import com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.peptide.Polypeptide;
import com.ariweiland.biophysics.peptide.Residue;

import java.util.Queue;

/**
 * @author Ari Weiland
 */
public abstract class Modeler {

    public static final int MAX_HEAP_SIZE = 4194304; // 262144, 524288, 1048576, 2097152, 4194304

    private final int dimension;

    protected Modeler(int dimension) {
        if (dimension < 2 || dimension > 3) {
            throw new IllegalArgumentException("Dimension of less than 2 or more than 3 does not make sense");
        }
        this.dimension = dimension;
    }

    public int getDimension() {
        return dimension;
    }

    /**
     * This method should stop the current folding process, if one is occurring.
     */
    public abstract void terminate();

    /**
     * This method takes in a polypeptide and returns a folded lattice.
     * @return
     */
    public abstract Lattice fold(Polypeptide polypeptide);

    /**
     * The core of the algorithm. This method should pop one Folding from the
     * queue. If that Folding is a solution state, it returns it, otherwise,
     * it should produce new Foldings and add them to the queue.
     * @param polypeptide
     * @param queue
     * @return
     */
    public abstract Folding iterate(Polypeptide polypeptide, Queue<Folding> queue);

    /**
     * For 2 dimensions, returns the perimeter of the
     * smallest rectangle the polypeptide can fit in.
     * For m^2 < n <= (m+1)^2,
     *  - returns 4m + 6 if n <= m(m+1)
     *  - returns 4m + 8 otherwise
     *
     * For 3 dimensions, returns the surface area of
     * the smallest box the polypeptide can fit in.
     * It returns 6 * m^2 where m = n^(1/3) + 1
     *
     * @param polypeptide
     * @return
     */
    protected int getSurfaceBound(Polypeptide polypeptide) {
        int n = polypeptide.size();
        if (getDimension() == 2) {
            int m = (int) Math.sqrt(n - 1);
            int maxPerim = 4 * m + 2;
            if (n > m * (m + 1)) {
                maxPerim += 2;
            }
            // add 4 because the ideal perimeter bound is overly limiting
            return maxPerim + 4;
        } else {
            // add 1 because the ideal surface bound is overly limiting
            double m = Math.pow(n, 1.0 / 3.0) + 1;
            return (int) (6 * m * m);
        }
    }

    /**
     * In the iteration phase, we can only account for favorable water interactions, because removing them
     * increases the energy value. Unfavorable water interactions, when removed, would decrease the energy
     * value, and that would cause the minimum to be an inaccurate heuristic.
     * @param p
     * @return
     */
    protected double getFavorableWaterInteraction(Peptide p) {
        return Math.min(p.interaction(Residue.H2O), 0);
    }

    public static void main(String[] args) {
        Modeler modeler = new CurrentParallelModeler(3);
//        Polypeptide polypeptide = Polypeptide.GLUCAGON;
//        Polypeptide polypeptide = new Polypeptide("(H)-(P)-(P)-(P)-(P)-(H)-(P)-(H)-(H)-(P)-(H)-(P)");
        Polypeptide polypeptide = new Polypeptide();
        for (int i=0; i<20; i++) {
            if (RandomUtils.tryChance(0.4)) {
                polypeptide.add(Residue.H);
            } else {
                polypeptide.add(Residue.P);
            }
        }
        System.out.println(polypeptide);
        System.out.println("Node count: " + polypeptide.size());
        System.out.println("Perimeter Bound: " + modeler.getSurfaceBound(polypeptide));
        System.out.println();

        long start = System.currentTimeMillis();
        Lattice lattice = modeler.fold(polypeptide);
        long elapsed = System.currentTimeMillis() - start;
        lattice.visualize();
        System.out.println("Elapsed time: " + (elapsed / 1000.0) + " s");
        System.out.println("Lattice energy: " + lattice.getEnergy());
        System.out.println("Perimeter: " + lattice.getSurfaceSize() + "/" + lattice.boundingPerimeter());
    }
}
