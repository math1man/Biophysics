package com.ariweiland.biophysics;

import java.util.concurrent.PriorityBlockingQueue;

/**
* @author Ari Weiland
*/
class PThread extends Thread {

    private final Polypeptide polypeptide;
    private final PriorityBlockingQueue<PState> initialHeap;
    private final PriorityBlockingQueue<PState> solutions;
    private final FixedHeap<PState> heap;

    PThread(Polypeptide polypeptide, PriorityBlockingQueue<PState> initialHeap,
            PriorityBlockingQueue<PState> solutions, int heapSize) {
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
            PState state = Modeler.iterate(polypeptide, heap);
            if (state != null) { // found a solution
                // don't bother with the solution if it isn't better than the current best
                // this will help conserve memory for larger polypeptides
                if (isWorthExploring(state)) {
                    solutions.put(state);
                }
                heap.clear();
            }
            if (heap.isEmpty()) {
                PState next = initialHeap.poll();
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

    public boolean isWorthExploring(PState next) {
        return solutions.isEmpty() || next.compareTo(solutions.peek()) < 0;
    }
}
