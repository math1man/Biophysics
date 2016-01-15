package com.ariweiland.biophysics.modeler;

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
     * Returns the perimeter of the smallest
     * rectangle with an area of at least n.
     * For m^2 < n <= (m+1)^2,
     *  - returns 4m + 6 if n <= m(m+1)
     *  - returns 4m + 8 otherwise
     * @param polypeptide
     * @return
     */
    protected int getPerimeterBound(Polypeptide polypeptide) {
        int n = polypeptide.size();
        int m = (int) Math.sqrt(n-1);
        int maxPerim = 4 * m + 2;
        if (n > m * (m+1)) {
            maxPerim += 2;
        }
        // add 4 because the ideal perimeter bound is overly limiting
        return maxPerim + 4;
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
        Modeler modeler = new CurrentSurfaceModeler(Residue.P);
//        Polypeptide polypeptide = Polypeptide.GLUCAGON;
//        Polypeptide polypeptide = new Polypeptide("(H)-(P)-(P)-(P)-(P)-(H)-(P)-(H)-(H)-(P)-(H)-(P)");
        Polypeptide polypeptide = new Polypeptide();
        for (int i=0; i<12; i++) {
            if (Math.random() < 0.4) {
                polypeptide.add(Residue.H);
            } else {
                polypeptide.add(Residue.P);
            }
        }
        System.out.println(polypeptide);
        System.out.println("Node count: " + polypeptide.size());
        System.out.println("Perimeter Bound: " + modeler.getPerimeterBound(polypeptide));
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
