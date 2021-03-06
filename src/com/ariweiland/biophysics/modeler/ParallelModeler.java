package com.ariweiland.biophysics.modeler;

import com.ariweiland.biophysics.lattice.Folding;
import com.ariweiland.biophysics.lattice.CheckedLattice;
import com.ariweiland.biophysics.peptide.Polypeptide;

import java.util.concurrent.PriorityBlockingQueue;
import java.util.concurrent.atomic.AtomicBoolean;

/**
 * @author Ari Weiland
 */
public abstract class ParallelModeler extends Modeler {

    private AtomicBoolean running = new AtomicBoolean();
    private PeptideThread[] threads;

    protected ParallelModeler(int dimension) {
        super(dimension);
    }

    /**
     * Dynamically calculates an ideal seed count for a given polypeptide.
     *
     * It uses the formula 10^(polypeptide.size() + 2).
     *
     * This results in a seed count of 1000 for size 10, 10000 for size 20,
     * 100000 for size 30, etc. with intermediate sizes being in between.
     *
     * @param polypeptide
     * @return
     */
    protected int getSeedCount(Polypeptide polypeptide) {
        double exponent = polypeptide.size() / 10.0 + 1.0; // TODO + getDimension() - 1.0;
        return (int) Math.pow(10.0, exponent);
    }

    /**
     * This helper method should initialize the heap in such a way that it contains all
     * symmetrically unique initial foldings. From these foldings, any other derived
     * folding should be completely unique.
     * @param polypeptide
     * @return
     */
    protected abstract PriorityBlockingQueue<Folding> initializeHeap(Polypeptide polypeptide);

    @Override
    public void terminate() {
        if (running.getAndSet(false) && threads != null) {
            for (PeptideThread thread : threads) {
                if (thread != null) {
                    thread.terminate();
                }
            }
        }
    }

    @Override
    public CheckedLattice fold(Polypeptide polypeptide) {
        running.set(true);
        PriorityBlockingQueue<Folding> initialHeap = initializeHeap(polypeptide);

        // iterate a few times to make the initial heap bigger
        int count = 0;
        while (running.get() && count < getSeedCount(polypeptide)) {
            Folding solution = iterate(polypeptide, initialHeap);
            if (solution != null) {
                return solution.lattice;
            }
            count++;
        }

        int processors = Runtime.getRuntime().availableProcessors();
        threads = new PeptideThread[processors];
        PriorityBlockingQueue<Folding> solutions = new PriorityBlockingQueue<>(initialHeap.size());

        System.out.println("Processors: " + processors);
        System.out.println("Initial Heap Size: " + initialHeap.size());
        for (int i=0; i< processors; i++) {
            threads[i] = new PeptideThread(this, polypeptide, initialHeap, solutions, MAX_HEAP_SIZE / processors - 1);
            if (running.get()) {
                threads[i].start();
            }
        }
        for (int i=0; i< processors; i++) {
            try {
                threads[i].join();
                count += threads[i].getCount();
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }
        }
        System.out.println("States Visited: " + count);
        if (running.get()) {
            return solutions.poll().lattice;
        } else {
            return new CheckedLattice(getDimension());
        }
    }
}
