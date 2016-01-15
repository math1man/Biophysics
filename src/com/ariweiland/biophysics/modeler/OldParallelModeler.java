package com.ariweiland.biophysics.modeler;

import com.ariweiland.biophysics.lattice.Direction;
import com.ariweiland.biophysics.lattice.Point2D;
import com.ariweiland.biophysics.lattice.Folding;
import com.ariweiland.biophysics.lattice.Lattice2D;
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
        Lattice2D line = new Lattice2D();
        line.put(0, 0, first);

        if (size > 1) {
            Peptide second = polypeptide.get(1);
            line.put(1, 0, second);

            // fill the queue initially.  this removes symmetrical solutions
            // if size == 2, the for loop will be ignored and none of this will matter
            double lowerBound = polypeptide.getMinEnergy() - 2 * first.minInteraction() - 2 * second.minInteraction();
            for (int i=2; i<size; i++) {
                Peptide next = polypeptide.get(i);
                lowerBound -= 2 * next.minInteraction();
                Lattice2D bend = new Lattice2D(line);
                bend.put(i - 1, 1, next);
                if (i == size - 1) {
                    lowerBound = bend.getEnergy();
                }
                initialHeap.add(new Folding(bend, i - 1, 1, i, lowerBound));
                line.put(i, 0, next);
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
            for (Direction d : Direction.values2D()) {
                Point2D next = folding.lastPoint.getAdjacent(d);
                if (!folding.lattice.containsPoint(next)) {
                    Lattice2D l = new Lattice2D(folding.lattice);
                    l.put(next, p);
                    // though limiting the protein to the smallest possible rectangle is
                    // overly limiting, empirically it seems that limiting it to a rectangle
                    // of perimeter 4 larger does not seem to restrict the solution at all
                    if (l.boundingPerimeter() <= getPerimeterBound(polypeptide)) {
                        double bound = folding.energyBound - 2 * p.minInteraction();
                        if (nextIndex < size - 1) {
                            for (Direction d1 : Direction.values2D()) {
                                if (d1 != d.getReverse()) {
                                    bound += p.interaction(l.get(next.getAdjacent(d)));
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
