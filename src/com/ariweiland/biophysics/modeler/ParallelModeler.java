package com.ariweiland.biophysics.modeler;

import com.ariweiland.biophysics.lattice.Folding;
import com.ariweiland.biophysics.lattice.Lattice;
import com.ariweiland.biophysics.peptide.Polypeptide;

import java.util.concurrent.PriorityBlockingQueue;

/**
 * @author Ari Weiland
 */
public abstract class ParallelModeler extends Modeler {

    protected abstract PriorityBlockingQueue<Folding> initializeHeap(Polypeptide polypeptide);

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
    protected int getSeedCount(Polypeptide polypeptide) {
        double exponent = polypeptide.size() / 10.0 + 1.0;
        return (int) Math.pow(10.0, exponent);
    }

    @Override
    public Lattice fold(Polypeptide polypeptide) {
        PriorityBlockingQueue<Folding> initialHeap = initializeHeap(polypeptide);

        // iterate a few times to make the initial heap bigger
        int count = 0;
        while (count < getSeedCount(polypeptide)) {
            Folding solution = iterate(polypeptide, initialHeap);
            if (solution != null) {
                return solution.lattice;
            }
            count++;
        }

        int processors = Runtime.getRuntime().availableProcessors();
        PeptideThread[] threads = new PeptideThread[processors];
        PriorityBlockingQueue<Folding> solutions = new PriorityBlockingQueue<>();

        System.out.println("Processors: " + processors);
        System.out.println("Initial Heap Size: " + initialHeap.size());
        for (int i=0; i< processors; i++) {
            threads[i] = new PeptideThread(this, polypeptide, initialHeap, solutions, MAX_HEAP_SIZE / processors - 1);
            threads[i].start();
        }
        for (int i=0; i< processors; i++) {
            try {
                threads[i].join();
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }
        }
        return solutions.poll().lattice;
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
}
