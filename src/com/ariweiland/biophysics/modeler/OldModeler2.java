package com.ariweiland.biophysics.modeler;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.lattice.Folding;
import com.ariweiland.biophysics.lattice.Lattice;
import com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.peptide.Polypeptide;
import com.ariweiland.biophysics.peptide.Residue;
import com.ariweiland.biophysics.FixedHeap;
import com.ariweiland.biophysics.Point;

import java.util.Queue;
import java.util.concurrent.atomic.AtomicBoolean;

/**
 * @author Ari Weiland
 */
public class OldModeler2 extends Modeler {

    private AtomicBoolean running = new AtomicBoolean();

    public OldModeler2() {
        super(2);
    }

    @Override
    public void terminate() {
        running.set(false);
    }

    @Override
    public Lattice fold(Polypeptide polypeptide) {
        running.set(true);
        // initialize the lattices
        int size = polypeptide.size();
        Peptide first = polypeptide.get(0);
        Lattice line = new Lattice(2);
        line.put(new Point(0, 0), first);
        if (size == 1) {
            return line;
        }
        Peptide second = polypeptide.get(1);
        line.put(new Point(1, 0), second);
        if (size == 2) {
            return line;
        }

        // fill the queue
        FixedHeap<Folding> heap = new FixedHeap<>(MAX_HEAP_SIZE - 1);
        double lowerBound = polypeptide.getMinEnergy(2) - 2 * first.minInteraction() - 2 * second.minInteraction();
        for (int i=2; i<size; i++) {
            Peptide next = polypeptide.get(i);
            lowerBound -= 2 * next.minInteraction();
            Lattice bend = new Lattice(line);
            Point point = new Point(i - 1, 1);
            bend.put(point, next);
            if (i == size - 1) {
                lowerBound = bend.getEnergy();
            }
            heap.add(new Folding(bend, point, i, lowerBound));
            line.put(new Point(i, 0), next);
        }
        heap.add(new Folding(line, new Point(size - 1, 0), size - 1, lowerBound));

        // begin the iteration
        Folding solution = null;
        int count = 0;
        while (running.get() && solution == null) {
            solution = iterate(polypeptide, heap);
            count++;
            if (count % 100000 == 0) {
                System.out.println(count + " states visited, " + heap.size() + " states in queue");
            }
        }
        System.out.println(count + " states visited, " + heap.size() + " states left in queue");
        if (solution == null) {
            return new Lattice(2);
        } else {
            return solution.lattice;
        }
    }

    @Override
    public Folding iterate(Polypeptide polypeptide, Queue<Folding> queue) {
        int size = polypeptide.size();
        Folding folding = queue.poll();
        int nextIndex = folding.index + 1;
        if (nextIndex < size) {
            Peptide p = polypeptide.get(nextIndex);
            for (Direction d : Direction.values(2)) {
                Point next = folding.lastPoint.getAdjacent(d);
                if (!folding.lattice.containsPoint(next)) {
                    Lattice l = new Lattice(folding.lattice);
                    l.put(next, p);
                    // though limiting the protein to the smallest possible rectangle is
                    // overly limiting, empirically it seems that limiting it to a rectangle
                    // of perimeter 4 larger does not seem to restrict the solution at all
                    if (l.boundingPerimeter() <= getSurfaceBound(polypeptide)) {
                        double bound = folding.energyBound - 2 * p.minInteraction()
                                - l.get(next.getAdjacent(d.getReverse())).interaction(Residue.H2O);
                        if (nextIndex < size - 1) {
                            for (Direction d1 : Direction.values(2)) {
                                if (d1 != d.getReverse()) {
                                    bound += p.interaction(l.get(next.getAdjacent(d1)));
                                }
                            }
                        } else {
                            bound = l.getEnergy();
                        }
                        queue.add(new Folding(l, next, nextIndex, bound));
                    }
                }
            }
            return null;
        } else {
            return folding;
        }
    }

}
