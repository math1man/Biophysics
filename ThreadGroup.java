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
    private int maxPerimFudgeFactor = 4;
    private boolean usePerimBound = true;
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

    /**
     * Returns the perimeter of the smallest
     * rectangle with an area of at least n.
     * For m^2 < n <= (m+1)^2,
     *  - returns 4m + 2 if n <= m(m+1)
     *  - returns 4m + 4 otherwise
     * Adds on the maxPerimFudgeFactor before returning
     * @param n
     * @return
     */
    private int getPerimBound(int n) {
        int m = (int) Math.sqrt(n-1);
        int maxPerim = 4 * m + 2;
        if (n > m * (m+1)) {
            maxPerim += 2;
        }
        return maxPerim + maxPerimFudgeFactor;
    }

    public int getMaxPerimFudgeFactor() {
        return maxPerimFudgeFactor;
    }

    public void setMaxPerimFudgeFactor(int maxPerimFudgeFactor) {
        this.maxPerimFudgeFactor = maxPerimFudgeFactor;
    }

    public boolean usePerimBound() {
        return usePerimBound;
    }

    public void setUsePerimBound(boolean usePerimBound) {
        this.usePerimBound = usePerimBound;
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
            int size = polypeptide.size();
            int count = 0;
            boolean running = !initialHeap.isEmpty();
            while (running) {
                if (heap.isEmpty()) {
                    if (!initialHeap.isEmpty()) {
                        PState next = initialHeap.poll();
                        if (solutions.isEmpty() || next.compareTo(solutions.peek()) < 0) {
                            heap.add(next);
                        } else {
                            // all other states will be worse
                            initialHeap.clear();
                            running = false;
                        }
                    } else {
                        running = false;
                    }
                }
                if (!heap.isEmpty()) {
                    PState state = heap.poll();
                    int index = state.index + 1;
                    if (index < size) {
                        Peptide p = polypeptide.get(index);
                        Point point = state.point;
                        double bound = state.energyBound - 2 * p.minInteraction();
                        for (Point.Direction d : Point.Direction.values()) {
                            Point next = point.getAdjacent(d);
                            if (!state.lattice.containsPoint(next)) {
                                Lattice l = new Lattice(state.lattice);
                                l.put(next, p);
                                // though limiting the protein to the smallest possible rectangle is
                                // overly limiting, empirically it seems that limiting it to a rectangle
                                // of perimeter 4 larger does not seem to restrict the solution at all
                                if (!usePerimBound || l.getMaxPerim() <= getPerimBound(size)) {
                                    double lb;
                                    if (index < size - 1) {
                                        lb = bound - l.get(next.getAdjacent(d.getReverse())).interaction(null);
                                        if (l.containsPoint(next.getAdjacent(d))) {
                                            lb += p.interaction(l.get(next.getAdjacent(d)));
                                        }
                                        if (l.containsPoint(next.getAdjacent(d.getLeft()))) {
                                            lb += p.interaction(l.get(next.getAdjacent(d.getLeft())));
                                        }
                                        if (l.containsPoint(next.getAdjacent(d.getRight()))) {
                                            lb += p.interaction(l.get(next.getAdjacent(d.getRight())));
                                        }
                                    } else {
                                        lb = l.getEnergy();
                                    }
                                    heap.add(new PState(l, next, index, lb));
                                }
                            }
                        }
                    } else {
                        solutions.put(state);
                        heap.clear();
                    }
                }
                count++;
                if (count % 100000 == 0) {
                    System.out.println(getName() + ": " + count + " states visited, " + heap.size() + " states in heap");
                }
            }
        }
    }
}
