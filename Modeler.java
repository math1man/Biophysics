package com.ariweiland.biophysics;

import java.util.PriorityQueue;
import java.util.Queue;
import java.util.concurrent.PriorityBlockingQueue;

/**
 * Possibly add to the Heuristic some function of the perimeter?
 * Confirmations with lower perimeter are better.
 *
 * @author Ari Weiland
 */
public class Modeler {

    public static void main(String[] args) {
        Polypeptide polypeptide = new Polypeptide();
        for (int i=0; i<30; i++) {
            if (Math.random() < 0.4) {
                polypeptide.add(PType.H);
            } else {
                polypeptide.add(PType.P);
            }
        }
        System.out.println(polypeptide);
        System.out.println("Node count: " + polypeptide.size());
        System.out.println("Perimeter Bound: " + getPerimBound(polypeptide.size()));
        System.out.println();

        long start = System.currentTimeMillis();
        Lattice lattice = fold(polypeptide, 500000);
        long elapsed = System.currentTimeMillis() - start;
        lattice.visualize();
        System.out.println("Elapsed time: " + (elapsed / 1000.0) + " s");
        System.out.println("Lattice energy: " + lattice.getEnergy());
        System.out.println("Perimeter: " + lattice.getPerimeter() + "/" + lattice.getMaxPerim());
//        System.out.println();
//        start = System.currentTimeMillis();
//        lattice = modeler.fold(500);
//        elapsed = System.currentTimeMillis() - start;
//        lattice.visualize();
//        System.out.println("Elapsed time: " + (elapsed/1000.0) + " s");
//        System.out.println("Lattice energy: " + lattice.getEnergy());
//        System.out.println("Perimeter: " + lattice.getPerimeter() + "/" + lattice.getMaxPerim());
    }

    private static final int MAX_HEAP_SIZE = 4194304; // 262144, 524288, 1048576, 2097152, 4194304

    /**
     * Returns the perimeter of the smallest
     * rectangle with an area of at least n.
     * For m^2 < n <= (m+1)^2,
     *  - returns 4m + 6 if n <= m(m+1)
     *  - returns 4m + 8 otherwise
     * @param n
     * @return
     */
    private static int getPerimBound(int n) {
        int m = (int) Math.sqrt(n-1);
        int maxPerim = 4 * m + 2;
        if (n > m * (m+1)) {
            maxPerim += 2;
        }
        return maxPerim + 4;
    }

    @Deprecated
    public static Lattice foldOld1(Polypeptide polypeptide) {
        // initialize the lattices
        int size = polypeptide.size();
        Peptide first = polypeptide.get(0);
        Lattice line = new Lattice();
        line.put(0, 0, first);
        if (size == 1) {
            return line;
        }
        Peptide second = polypeptide.get(1);
        line.put(1, 0, second);
        if (size == 2) {
            return line;
        }

        // fill the queue
        PriorityQueue<PState> pq = new PriorityQueue<PState>();
        double lowerBound = polypeptide.getMinEnergy() - 2 * first.minInteraction() - 2 * second.minInteraction();
        for (int i=2; i<size; i++) {
            Peptide next = polypeptide.get(i);
            lowerBound -= 2 * next.minInteraction();
            Lattice bend = new Lattice(line);
            bend.put(i-1, 1, next);
            if (i == size-1) {
                lowerBound = bend.getEnergy();
            }
            pq.add(new PState(bend, i - 1, 1, i, lowerBound));
            line.put(i, 0, next);
        }
        pq.add(new PState(line, size - 1, 0, size - 1, lowerBound));

        // begin the iteration
        Lattice solution = null;
        int count = 0;
        while (solution == null) {
            PState state = pq.poll();
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
                        if (l.getMaxPerim() <= getPerimBound(size)) {
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
                            pq.add(new PState(l, next, index, lb));
                        }
                    }
                }
            } else {
                solution = state.lattice;
            }
            count++;
            if (count % 10000 == 0) {
                System.out.println(count + " states visited, " + pq.size() + " states in queue");
            }
        }
        System.out.println(count + " states visited, " + pq.size() + " states left in queue");
        return solution;
    }

    @Deprecated
    public static Lattice foldOld2(Polypeptide polypeptide) {
        // initialize the lattices
        int size = polypeptide.size();
        Peptide first = polypeptide.get(0);
        Lattice line = new Lattice();
        line.put(0, 0, first);
        if (size == 1) {
            return line;
        }
        Peptide second = polypeptide.get(1);
        line.put(1, 0, second);
        if (size == 2) {
            return line;
        }

        // fill the queue
        FixedHeap<PState> heap = new FixedHeap<PState>(MAX_HEAP_SIZE - 1);
        double lowerBound = polypeptide.getMinEnergy() - 2 * first.minInteraction() - 2 * second.minInteraction();
        for (int i=2; i<size; i++) {
            Peptide next = polypeptide.get(i);
            lowerBound -= 2 * next.minInteraction();
            Lattice bend = new Lattice(line);
            bend.put(i-1, 1, next);
            if (i == size-1) {
                lowerBound = bend.getEnergy();
            }
            heap.add(new PState(bend, i - 1, 1, i, lowerBound));
            line.put(i, 0, next);
        }
        heap.add(new PState(line, size - 1, 0, size - 1, lowerBound));

        // begin the iteration
        PState solution = null;
        int count = 0;
        while (solution == null) {
            solution = iterate(polypeptide, heap);
            count++;
            if (count % 100000 == 0) {
                System.out.println(count + " states visited, " + heap.size() + " states in queue");
            }
        }
        System.out.println(count + " states visited, " + heap.size() + " states left in queue");
        return solution.lattice;
    }

    /**
     * The fastest version of the algorithm yet. Parallelization works!
     * Useful for all sizes of Polypeptides, assuming an appropriate seed
     * count is used. Seed count should be kept small (~1000) but may need
     * to be increased so that the heaps don't fill up.
     *     *
     * @param seedCount
     * @return
     */
    public static Lattice fold(Polypeptide polypeptide, int seedCount) {
        // initialize the lattices
        int size = polypeptide.size();
        Peptide first = polypeptide.get(0);
        Lattice line = new Lattice();
        line.put(0, 0, first);
        if (size == 1) {
            return line;
        }
        Peptide second = polypeptide.get(1);
        line.put(1, 0, second);
        if (size == 2) {
            return line;
        }

        // fill the queue initially.  this removes symmetrical solutions
        PriorityBlockingQueue<PState> initialHeap = new PriorityBlockingQueue<PState>();
        double lowerBound = polypeptide.getMinEnergy() - 2 * first.minInteraction() - 2 * second.minInteraction();
        for (int i=2; i<size; i++) {
            Peptide next = polypeptide.get(i);
            lowerBound -= 2 * next.minInteraction();
            Lattice bend = new Lattice(line);
            bend.put(i-1, 1, next);
            if (i == size-1) {
                lowerBound = bend.getEnergy();
            }
            initialHeap.add(new PState(bend, i - 1, 1, i, lowerBound));
            line.put(i, 0, next);
        }
        initialHeap.add(new PState(line, size - 1, 0, size - 1, lowerBound));

        // iterate a few times to make the initial heap bigger
        int count = 0;
        while (count < seedCount) {
            PState solution = iterate(polypeptide, initialHeap);
            if (solution != null) {
                return solution.lattice;
            }
            count++;
        }

        ThreadGroup group = new ThreadGroup(polypeptide, initialHeap);
        group.setTotalHeapSize(MAX_HEAP_SIZE);
        return group.process();
    }

    public static PState iterate(Polypeptide polypeptide, Queue<PState> queue) {
        int size = polypeptide.size();
        PState state = queue.poll();
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
                    if (l.getMaxPerim() <= getPerimBound(size)) {
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
                        queue.add(new PState(l, next, index, lb));
                    }
                }
            }
            return null;
        } else {
            return state;
        }
    }
}
