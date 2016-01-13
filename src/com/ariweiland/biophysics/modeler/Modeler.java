package com.ariweiland.biophysics.src.com.ariweiland.biophysics.modeler;

import com.ariweiland.biophysics.src.com.ariweiland.biophysics.lattice.Folding;
import com.ariweiland.biophysics.src.com.ariweiland.biophysics.lattice.Lattice;
import com.ariweiland.biophysics.src.com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.src.com.ariweiland.biophysics.peptide.Polypeptide;
import com.ariweiland.biophysics.src.com.ariweiland.biophysics.peptide.Residue;

import java.util.Queue;
import java.util.concurrent.PriorityBlockingQueue;

/**
 * @author Ari Weiland
 */
public abstract class Modeler {

    public static final int MAX_HEAP_SIZE = 4194304; // 262144, 524288, 1048576, 2097152, 4194304

    /**
     * Seed count is used for parallelization. It should be kept small (~1000) but may need
     * to be increased so that the heaps don't fill up. The default constructor sets it to 0.
     */
    private final int seedCount;

    protected Modeler() {
        this(0);
    }

    protected Modeler(int seedCount) {
        this.seedCount = seedCount;
    }

    public int getSeedCount() {
        return seedCount;
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
        Modeler modeler = new SurfaceModeler(10000, Residue.POS);
//        Polypeptide polypeptide = new Polypeptide("+PP PHPP-HP+HH-P++HP-HHPHHHPP");
//        Polypeptide polypeptide = new Polypeptide("PHPHPPHPHPPP");
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
        Lattice lattice = modeler.fold(polypeptide);
        long elapsed = System.currentTimeMillis() - start;
        lattice.visualize();
        System.out.println("Elapsed time: " + (elapsed / 1000.0) + " s");
        System.out.println("Lattice energy: " + lattice.getEnergy());
        System.out.println("Perimeter: " + lattice.getPerimeter() + "/" + lattice.boundingPerimeter());
    }
}
