package com.ariweiland.biophysics;

import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.PriorityBlockingQueue;

/**
 * @author Ari Weiland
 */
public class ThreadGroup {

    /**
     * Though limiting the protein to the smallest possible
     * rectangle is overly restrictive, empirically it seems
     * that limiting it to a rectangle of perimeter 4 larger
     * does not seem to restrict the solution at all.
     */
    private int totalHeapSize = 4194304;

    private final Polypeptide polypeptide;
    private final PriorityBlockingQueue<PState> initialHeap;
    private final PriorityBlockingQueue<PState> solutions = new PriorityBlockingQueue<PState>();

    public ThreadGroup(Polypeptide polypeptide, PriorityBlockingQueue<PState> initialHeap) {
        this.polypeptide = polypeptide;
        this.initialHeap = initialHeap;
    }

    public Lattice process() {
        int processors = Runtime.getRuntime().availableProcessors();
        System.out.println("Processors: " + processors);
        int heapSize = totalHeapSize / processors - 1;
        List<PThread> threads = new LinkedList<PThread>();
        for (int i=0; i<processors; i++) {
            PThread thread = new PThread(heapSize);
            threads.add(thread);
            thread.start();
        }
        for (int i=0; i<processors; i++) {
            try {
                threads.get(i).join();
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }
        }
        return solutions.poll().lattice;
    }

    public int getTotalHeapSize() {
        return totalHeapSize;
    }

    public void setTotalHeapSize(int totalHeapSize) {
        this.totalHeapSize = totalHeapSize;
    }

    private class PThread extends Thread {

        private final FixedHeap<PState> heap;

        private PThread(int heapSize) {
            this.heap = new FixedHeap<PState>(heapSize);
        }

        @Override
        public void run() {
            super.run();
            int count = 0;
            boolean running = !initialHeap.isEmpty();
            while (running) {
                if (heap.isEmpty()) {
                    PState next = initialHeap.poll();
                    if (next != null && (solutions.isEmpty() || next.compareTo(solutions.peek()) < 0)) {
                        heap.add(next);
                    } else {
                        // all other states will be worse
                        initialHeap.clear();
                        running = false;
                    }
                }
                if (!heap.isEmpty()) {
                    PState state = Modeler.iterate(polypeptide, heap);
                    if (state != null) {
                        // don't bother with the solution if it isn't better than the current best
                        // this will help conserve memory for larger polypeptides
                        if (solutions.isEmpty() || state.compareTo(solutions.peek()) < 0) {
                            solutions.put(state);
                        }
                        heap.clear();
                    }
                }
                count++;
                if (count % 1000000 == 0) {
                    System.out.println(getName() + ": " + count/1000000 + "M states visited, " + heap.size() + " states in heap, " + initialHeap.size() + " states left in initial heap");
                }
            }
        }
    }
}
