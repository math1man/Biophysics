package com.ariweiland.biophysics.modeler;

import com.ariweiland.biophysics.lattice.Folding;
import com.ariweiland.biophysics.lattice.Lattice;
import com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.peptide.Polypeptide;
import com.ariweiland.biophysics.peptide.Residue;

import java.util.Queue;
import java.util.concurrent.PriorityBlockingQueue;

/**
 * @author Ari Weiland
 */
public abstract class Modeler {

    public static final int MAX_HEAP_SIZE = 4194304; // 262144, 524288, 1048576, 2097152, 4194304

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
     * Primary method that generates a folded protein from the
     * @param processors
     * @return
     */
    public Lattice parallelize(Polypeptide polypeptide, PriorityBlockingQueue<Folding> initialHeap,
                               int processors, int processHeapSize) {
        System.out.println("Processors: " + processors);
        System.out.println("Initial Heap Size: " + initialHeap.size());
        PeptideThread[] threads = new PeptideThread[processors];
        PriorityBlockingQueue<Folding> solutions = new PriorityBlockingQueue<>();

        for (int i=0; i<processors; i++) {
            threads[i] = new PeptideThread(this, polypeptide, initialHeap, solutions, processHeapSize);
            threads[i].start();
        }
        for (int i=0; i<processors; i++) {
            try {
                threads[i].join();
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }
        }
        return solutions.poll().lattice;
    }

    /**
     * Dynamically calculates an ideal seed count for a given polypeptide.
     *
     * It uses the formula 10^(polypeptide.size() + 2).
     *
     * This results in a seed count of 1000 for size 10, 10000 for size 20,
     * 100000 for size 30, etc. with intermediate sizes being in between.
     * @param polypeptide
     * @return
     */
    public static int getSeedCount(Polypeptide polypeptide) {
        double exponent = polypeptide.size() / 10.0 + 2.0;
        return (int) Math.pow(10.0, exponent);
    }

    /**
     * Returns the perimeter of the smallest
     * rectangle with an area of at least n.
     * For m^2 < n <= (m+1)^2,
     *  - returns 4m + 6 if n <= m(m+1)
     *  - returns 4m + 8 otherwise
     * @param polypeptide
     * @return
     */
    public static int getPerimeterBound(Polypeptide polypeptide) {
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
     * This method calculates the maximum y-value a polypeptide should ever reach in surface modeling.
     * It is related to the perimeter bound.
     * @param polypeptide
     * @return
     */
    public static int getMaxY(Polypeptide polypeptide) {
        return getPerimeterBound(polypeptide) / 4 + 2;
    }

    /**
     * In the iteration phase, we can only account for favorable water interactions, because removing them
     * increases the energy value. Unfavorable water interactions, when removed, would decrease the energy
     * value, and that would cause the minimum to be an inaccurate heuristic.
     * @param p
     * @return
     */
    public static double getFavorableWaterInteraction(Peptide p) {
        return Math.min(p.interaction(Residue.H2O), 0);
    }

    public static void main(String[] args) {
        Modeler modeler = new SurfaceModeler(Residue.P);
//        Polypeptide polypeptide = Polypeptide.GLUCAGON;
//        Polypeptide polypeptide = new Polypeptide("(H)-(P)-(P)-(P)-(P)-(H)-(P)-(H)-(H)-(P)-(H)-(P)");
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
        System.out.println("Perimeter Bound: " + getPerimeterBound(polypeptide));
        System.out.println();

        long start = System.currentTimeMillis();
        Lattice lattice = modeler.fold(polypeptide);
        long elapsed = System.currentTimeMillis() - start;
        lattice.visualize();
        System.out.println("Elapsed time: " + (elapsed / 1000.0) + " s");
        System.out.println("Lattice energy: " + lattice.getEnergy());
        System.out.println("Perimeter: " + lattice.getPerimeter() + "/" + lattice.boundingPerimeter());
    }
}
