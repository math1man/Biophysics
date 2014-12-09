package com.ariweiland.biophysics;

import java.util.ArrayList;
import java.util.List;
import java.util.PriorityQueue;
import java.util.concurrent.PriorityBlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * @author Ari Weiland
 */
public class Modeler {

    public static Lattice fold(Polypeptide polypeptide) {
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
                double bound = state.lowerBound - 2 * p.minInteraction();
                for (Point.Direction d : Point.Direction.values()) {
                    Point next = point.getAdjacent(d);
                    if (!state.lattice.containsPoint(next)) {
                        Lattice l = new Lattice(state.lattice);
                        l.put(next, p);
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
            } else {
                solution = state.lattice;
            }
            count++;
            if (count % 10000 == 0) {
                System.out.println(count + " states visited, " + pq.size() + " states in queue");
            }
        }
        return solution;
    }

    public static Lattice foldParallelized(final Polypeptide polypeptide) {
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
                        double bound = state.lowerBound - 2 * p.minInteraction();
                        for (Point.Direction d : Point.Direction.values()) {
                            Point next = point.getAdjacent(d);
                            if (!state.lattice.containsPoint(next)) {
                                Lattice l = new Lattice(state.lattice);
                                l.put(next, p);
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
                                pbq.put(new Modeler.State(l, next, index, lb));
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
        final double lowerBound;

        private State(Lattice lattice, int x, int y, int index, double lowerBound) {
            this(lattice, new Point(x, y), index, lowerBound);
        }

        private State(Lattice lattice, Point point, int index, double lowerBound) {
            this.lattice = lattice;
            this.point = point;
            this.index = index;
            this.lowerBound = lowerBound;
        }

        @Override
        public int compareTo(State o) {
            return Double.compare(lowerBound, o.lowerBound);
        }

        @Override
        public String toString() {
            return "" + lowerBound;
        }
    }

    public static void main(String[] args) {
        Polypeptide polypeptide = new Polypeptide();
        for (int i=0; i<18; i++) {
            if (Math.random() < 0.4) {
                polypeptide.add(Type.NONPOLAR);
            } else {
                polypeptide.add(Type.POLAR);
            }
        }
        System.out.println(polypeptide);
        long start = System.currentTimeMillis();
        Lattice lattice = foldParallelized(polypeptide);
        long elapsed = System.currentTimeMillis() - start;
        lattice.visualize();
        System.out.println("Elapsed time: " + (elapsed/1000.0) + " s");
        System.out.println("Lattice energy: " + lattice.getEnergy());
    }
}
