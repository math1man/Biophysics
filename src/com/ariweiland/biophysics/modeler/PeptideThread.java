package com.ariweiland.biophysics.src.com.ariweiland.biophysics.modeler;

import com.ariweiland.biophysics.src.com.ariweiland.biophysics.FixedHeap;
import com.ariweiland.biophysics.src.com.ariweiland.biophysics.lattice.Folding;
import com.ariweiland.biophysics.src.com.ariweiland.biophysics.peptide.Polypeptide;

import java.util.concurrent.PriorityBlockingQueue;

/**
 * This thread implementation is used to parallelize Modeler's iterate method.
 * It builds its own heap starting with Foldings it pulls from the initialHeap,
 * and if it finds a solution, adds it to solutions and pulls a new seed from
 * the initialHeap to try.
 * @author Ari Weiland
 */
public class PeptideThread extends Thread {

    private final Modeler modeler;
    private final Polypeptide polypeptide;
    private final PriorityBlockingQueue<Folding> initialHeap;
    private final PriorityBlockingQueue<Folding> solutions;
    private final FixedHeap<Folding> heap;

    public PeptideThread(Modeler modeler, Polypeptide polypeptide, PriorityBlockingQueue<Folding> initialHeap,
                         PriorityBlockingQueue<Folding> solutions, int heapSize) {
        this.modeler = modeler;
        this.polypeptide = polypeptide;
        this.initialHeap = initialHeap;
        this.solutions = solutions;
        this.heap = new FixedHeap<>(heapSize);
    }

    @Override
    public void run() {
        int count = 0;
        if (!initialHeap.isEmpty()) {
            heap.add(initialHeap.poll());
        }
        while (!heap.isEmpty()) {
            Folding state = modeler.iterate(polypeptide, heap);
            if (state != null) { // found a solution
                // don't bother with the solution if it isn't better than the current best
                // this will help conserve memory for larger polypeptides
                if (isWorthExploring(state)) {
                    solutions.put(state);
                }
                heap.clear();
            }
            if (heap.isEmpty()) {
                Folding next = initialHeap.poll();
                if (next != null && isWorthExploring(next)) {
                    heap.add(next);
                }
            }
            count++;
            if (count % 1000000 == 0) {
                System.out.println(getName() + ": " + count/1000000 + "M states visited, " + heap.size() + " states in heap, " + initialHeap.size() + " states left in initial heap");
            }
        }
    }

    public boolean isWorthExploring(Folding next) {
        return solutions.isEmpty() || next.compareTo(solutions.peek()) < 0;
    }
}
