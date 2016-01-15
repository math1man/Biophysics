package com.ariweiland.biophysics.modeler;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.lattice.Folding;
import com.ariweiland.biophysics.lattice.Lattice;
import com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.peptide.Polypeptide;

import java.util.Queue;
import java.util.concurrent.PriorityBlockingQueue;

/**
 * Possibly add to the Heuristic some function of the perimeter?
 * Confirmations with lower perimeter are better.
 *
 * @author Ari Weiland
 */
public class OldParallelModeler extends ParallelModeler {

    @Override
    protected PriorityBlockingQueue<Folding> initializeHeap(Polypeptide polypeptide) {
        PriorityBlockingQueue<Folding> initialHeap = new PriorityBlockingQueue<>();
        int size = polypeptide.size();
        // initialize the lattices
        Peptide first = polypeptide.get(0);
        Lattice line = new Lattice(2);
        line.put(first, 0, 0);

        if (size > 1) {
            Peptide second = polypeptide.get(1);
            line.put(second, 1, 0);

            // fill the queue initially.  this removes symmetrical solutions
            // if size == 2, the for loop will be ignored and none of this will matter
            double lowerBound = polypeptide.getMinEnergy(2) - 2 * first.minInteraction() - 2 * second.minInteraction();
            for (int i=2; i<size; i++) {
                Peptide next = polypeptide.get(i);
                lowerBound -= 2 * next.minInteraction();
                Lattice bend = new Lattice(line);
                bend.put(next, i - 1, 1);
                if (i == size - 1) {
                    lowerBound = bend.getEnergy();
                }
                initialHeap.add(new Folding(bend, i - 1, 1, i, lowerBound));
                line.put(next, i, 0);
            }
        }
        initialHeap.add(new Folding(line, size - 1, 0, size - 1, line.getEnergy()));
        return initialHeap;
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
                    l.put(p, next);
                    // though limiting the protein to the smallest possible rectangle is
                    // overly limiting, empirically it seems that limiting it to a rectangle
                    // of perimeter 4 larger does not seem to restrict the solution at all
                    if (l.boundingPerimeter() <= getPerimeterBound(polypeptide)) {
                        double bound = folding.energyBound - 2 * p.minInteraction();
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
