package com.ariweiland.biophysics;

import java.util.ArrayList;
import java.util.List;
import java.util.PriorityQueue;
import java.util.concurrent.PriorityBlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Possibly add to the Heuristic some function of the perimeter?
 * Confirmations with lower perimeter are better.
 *
 * @author Ari Weiland
 */
public class Modeler {

    /**
     * Though limiting the protein to the smallest possible
     * rectangle is overly restrictive, empirically it seems
     * that limiting it to a rectangle of perimeter 4 larger
     * does not seem to restrict the solution at all.
     */
    public static int MAX_PERIM_FUDGE_FACTOR = 4;
    public static boolean USE_MAX_PERIM = true;
    public static boolean USE_PERIMETER = true;

    /**
     * Returns the perimeter of the smallest
     * rectangle with an area of at least n.
     * For m^2 < n <= (m+1)^2,
     *  - returns 4m + 2 if n <= m(m+1)
     *  - returns 4m + 4 otherwise
     * Adds on the MAX_PERIM_FUDGE_FACTOR before returning
     * @param n
     * @return
     */
    private static int getPerimBound(int n) {
        int m = (int) Math.sqrt(n-1);
        int maxPerim = 4 * m + 2;
        if (n > m * (m+1)) {
            maxPerim += 2;
        }
        return maxPerim + MAX_PERIM_FUDGE_FACTOR;
    }

    public static Lattice fold(Polypeptide polypeptide) {
        // the perimeter preference is generally unproductive
        USE_PERIMETER = false;
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
        PriorityQueue<State> pq = new PriorityQueue<State>(6000000);
        double lowerBound = polypeptide.getMinEnergy() - 2 * first.minInteraction() - 2 * second.minInteraction();
        for (int i=2; i<size; i++) {
            Peptide next = polypeptide.get(i);
            lowerBound -= 2 * next.minInteraction();
            Lattice bend = new Lattice(line);
            bend.put(i-1, 1, next);
            if (i == size-1) {
                lowerBound = bend.getEnergy();
            }
            pq.add(new State(bend, i-1, 1, i, lowerBound));
            line.put(i, 0, next);
        }
        pq.add(new State(line, size-1, 0, size-1, lowerBound));

        // begin the iteration
        Lattice solution = null;
        int count = 0;
        while (solution == null) {
            State state = pq.poll();
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
                        if (!USE_MAX_PERIM || l.getMaxPerim() <= getPerimBound(size)) {
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
                            pq.add(new State(l, next, index, lb));
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

    public static Lattice foldLimited(Polypeptide polypeptide) {
        // need as much preference as we can get
        // since we are cutting off states
        USE_PERIMETER = true;
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
        FixedHeap<State> heap = new FixedHeap<State>(4194303); // 262143, 524287, 1048575, 2097151, 4194303
        double lowerBound = polypeptide.getMinEnergy() - 2 * first.minInteraction() - 2 * second.minInteraction();
        for (int i=2; i<size; i++) {
            Peptide next = polypeptide.get(i);
            lowerBound -= 2 * next.minInteraction();
            Lattice bend = new Lattice(line);
            bend.put(i-1, 1, next);
            if (i == size-1) {
                lowerBound = bend.getEnergy();
            }
            heap.push(new State(bend, i - 1, 1, i, lowerBound));
            line.put(i, 0, next);
        }
        heap.push(new State(line, size - 1, 0, size - 1, lowerBound));

        // begin the iteration
        Lattice solution = null;
        int count = 0;
        while (solution == null) {
            State state = heap.pop();
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
                        if (!USE_MAX_PERIM || l.getMaxPerim() <= getPerimBound(size)) {
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
                            heap.push(new State(l, next, index, lb));
                        }
                    }
                }
            } else {
                solution = state.lattice;
            }
            count++;
            if (count % 10000 == 0) {
                System.out.println(count + " states visited, " + heap.size() + " states in queue");
            }
        }
        System.out.println(count + " states visited, " + heap.size() + " states left in queue");
        return solution;
    }

    public static Lattice foldParallelized(final Polypeptide polypeptide) {
        // the perimeter preference is generally unproductive
        USE_PERIMETER = false;
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
        PriorityBlockingQueue<State> pbq = new PriorityBlockingQueue<State>(6000000);
        double lowerBound = polypeptide.getMinEnergy() - 2 * first.minInteraction() - 2 * second.minInteraction();
        for (int i=2; i<size; i++) {
            Peptide next = polypeptide.get(i);
            lowerBound -= 2 * next.minInteraction();
            Lattice bend = new Lattice(line);
            bend.put(i-1, 1, next);
            if (i == size-1) {
                lowerBound = bend.getEnergy();
            }
            pbq.put(new State(bend, i - 1, 1, i, lowerBound));
            line.put(i, 0, next);
        }
        pbq.put(new State(line, size - 1, 0, size - 1, lowerBound));

        // begin parallelization steps
        int processors = Runtime.getRuntime().availableProcessors();
        System.out.println(processors);
        List<PeptideThread> threads = new ArrayList<PeptideThread>();
        for (int i=0; i<processors; i++) {
            PeptideThread thread = new PeptideThread(polypeptide, pbq);
            threads.add(thread);
            thread.start();
        }

        // collect solutions
        List<Lattice> solutions = new ArrayList<Lattice>();
        for (int i=0; i<processors; i++) {
            try {
                threads.get(i).join();
                solutions.add(threads.get(i).getSolution());
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }
        }

        // find best solution
        Lattice solution = solutions.get(0);
        double energy = solution.getEnergy();
        for (Lattice l : solutions) {
            if (l.getEnergy() < energy) {
                solution = l;
                energy = l.getEnergy();
            }
        }
        return solution;
    }

    private static class PeptideThread extends Thread {

        public static AtomicInteger count = new AtomicInteger(0);

        private final Polypeptide polypeptide;
        private final PriorityBlockingQueue<Modeler.State> pbq;
        private Lattice solution = null;

        private PeptideThread(Polypeptide polypeptide, PriorityBlockingQueue<Modeler.State> pbq) {
            this.polypeptide = polypeptide;
            this.pbq = pbq;
        }

        public Lattice getSolution() {
            return solution;
        }

        @Override
        public void run() {
            super.run();
            int size = polypeptide.size();
            while (solution == null) {
                try {
                    Modeler.State state = pbq.take();
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
                                if (!USE_MAX_PERIM || l.getMaxPerim() <= getPerimBound(size)) {
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
                                    pbq.add(new Modeler.State(l, next, index, lb));
                                }
                            }
                        }
                    } else {
                        solution = state.lattice;
                        pbq.offer(state);
                    }
                    int temp = count.incrementAndGet();
                    if (temp % 10000 == 0) {
                        System.out.println(temp + " states visited, " + pbq.size() + " states in queue");
                    }
                } catch (InterruptedException e) {
                    throw new RuntimeException(e);
                }
            }
        }
    }

    private static class State implements Comparable<State> {
        final Lattice lattice;
        final Point point;
        final int index;
        final double energyBound;
        final int perimeter;

        private State(Lattice lattice, int x, int y, int index, double energyBound) {
            this(lattice, new Point(x, y), index, energyBound);
        }

        private State(Lattice lattice, Point point, int index, double energyBound) {
            this.lattice = lattice;
            this.point = point;
            this.index = index;
            this.energyBound = energyBound;
            this.perimeter = lattice.getPerimeter();
        }

        @Override
        public int compareTo(State o) {
            int compare = Double.compare(energyBound, o.energyBound);
            if (compare == 0 && USE_PERIMETER) {
                compare = Integer.compare(perimeter, o.perimeter);
            }
            return compare;
        }

        @Override
        public String toString() {
            return energyBound + "/" + perimeter;
        }
    }

    public static void main(String[] args) {
        Polypeptide polypeptide = new Polypeptide();
        for (int i=0; i<24; i++) {
            if (Math.random() < 0.4) {
                polypeptide.add(Type.H);
            } else {
                polypeptide.add(Type.P);
            }
        }
        USE_MAX_PERIM = true;
        System.out.println(polypeptide);
        System.out.println("Node count: " + polypeptide.size());
        System.out.println("Perimeter Bound: " + getPerimBound(polypeptide.size()));
        long start = System.currentTimeMillis();
        Lattice lattice = foldLimited(polypeptide);
        long elapsed = System.currentTimeMillis() - start;
        lattice.visualize();
        System.out.println("Elapsed time: " + (elapsed/1000.0) + " s");
        System.out.println("Lattice energy: " + lattice.getEnergy());
        System.out.println("Perimeter: " + lattice.getPerimeter() + "/" + lattice.getMaxPerim());
//        System.out.println();
//        USE_MAX_PERIM = false;
//        System.out.println("Not using max perimeter");
        start = System.currentTimeMillis();
        lattice = fold(polypeptide);
        elapsed = System.currentTimeMillis() - start;
        lattice.visualize();
        System.out.println("Elapsed time: " + (elapsed/1000.0) + " s");
        System.out.println("Lattice energy: " + lattice.getEnergy());
        System.out.println("Perimeter: " + lattice.getPerimeter() + "/" + lattice.getMaxPerim());
    }
}
