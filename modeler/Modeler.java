package com.ariweiland.biophysics.modeler;

import com.ariweiland.biophysics.lattice.Folding;
import com.ariweiland.biophysics.peptide.Polypeptide;
import com.ariweiland.biophysics.peptide.Residue;
import com.ariweiland.biophysics.lattice.Lattice;

import java.util.Queue;

/**
 * @author Ari Weiland
 */
public abstract class Modeler {

    public static final int MAX_HEAP_SIZE = 4194304; // 262144, 524288, 1048576, 2097152, 4194304

    /**
     * Returns the perimeter of the smallest
     * rectangle with an area of at least n.
     * For m^2 < n <= (m+1)^2,
     *  - returns 4m + 6 if n <= m(m+1)
     *  - returns 4m + 8 otherwise
     * @param n
     * @return
     */
    public static int getPerimBound(int n) {
        int m = (int) Math.sqrt(n-1);
        int maxPerim = 4 * m + 2;
        if (n > m * (m+1)) {
            maxPerim += 2;
        }
        // add 4 because the ideal perimeter bound is overly limiting
        return maxPerim + 4;
    }

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

    public static void main(String[] args) {
        Polypeptide polypeptide = new Polypeptide();
        for (int i=0; i<20; i++) {
            if (Math.random() < 0.4) {
                polypeptide.add(Residue.H);
            } else {
                polypeptide.add(Residue.P);
            }
        }
        System.out.println(polypeptide);
        System.out.println("Node count: " + polypeptide.size());
        System.out.println("Perimeter Bound: " + getPerimBound(polypeptide.size()));
        System.out.println();

        long start = System.currentTimeMillis();
        Lattice lattice = new SurfaceModeler(1000, Residue.H).fold(polypeptide);
        long elapsed = System.currentTimeMillis() - start;
        lattice.visualize();
        System.out.println("Elapsed time: " + (elapsed / 1000.0) + " s");
        System.out.println("Lattice energy: " + lattice.getEnergy());
        System.out.println("Perimeter: " + lattice.getPerimeter() + "/" + lattice.boundingPerimeter());
    }

}
