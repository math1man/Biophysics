package com.ariweiland.biophysics.modeler;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.lattice.Folding;
import com.ariweiland.biophysics.lattice.Lattice;
import com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.peptide.Polypeptide;
import com.ariweiland.biophysics.peptide.Residue;
import com.ariweiland.biophysics.Point;

import java.util.PriorityQueue;
import java.util.Queue;
import java.util.concurrent.atomic.AtomicBoolean;

/**
 * @author Ari Weiland
 */
public class OldModeler1 extends Modeler {

    private AtomicBoolean running = new AtomicBoolean();

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
        PriorityQueue<Folding> pq = new PriorityQueue<>();
        double lowerBound = polypeptide.getMinEnergy() - 2 * first.minInteraction() - 2 * second.minInteraction();
        for (int i=2; i<size; i++) {
            Peptide next = polypeptide.get(i);
            lowerBound -= 2 * next.minInteraction();
            Lattice bend = new Lattice(line);
            bend.put(i-1, 1, next);
            if (i == size-1) {
                lowerBound = bend.getEnergy();
            }
            pq.add(new Folding(bend, i - 1, 1, i, lowerBound));
            line.put(i, 0, next);
        }
        pq.add(new Folding(line, size - 1, 0, size - 1, lowerBound));

        // begin the iteration
        Folding solution = null;
        int count = 0;
        while (running.get() && solution == null) {
            solution = iterate(polypeptide, pq);
            count++;
            if (count % 10000 == 0) {
                System.out.println(count + " states visited, " + pq.size() + " states in queue");
            }
        }
        System.out.println(count + " states visited, " + pq.size() + " states left in queue");
        if (solution == null) {
            return new Lattice();
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
            for (Direction d : Direction.values()) {
                Point next = folding.lastPoint.getAdjacent(d);
                if (!folding.lattice.containsPoint(next)) {
                    Lattice l = new Lattice(folding.lattice);
                    l.put(next, p);
                    // though limiting the protein to the smallest possible rectangle is
                    // overly limiting, empirically it seems that limiting it to a rectangle
                    // of perimeter 4 larger does not seem to restrict the solution at all
                    if (l.boundingPerimeter() <= getPerimeterBound(polypeptide)) {
                        double bound = folding.energyBound - 2 * p.minInteraction()
                                - l.get(next.getAdjacent(d.getReverse())).interaction(Residue.H2O);
                        if (nextIndex < size - 1) {
                            for (Direction d1 : Direction.values()) {
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
