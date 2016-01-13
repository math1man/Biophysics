package com.ariweiland.biophysics.modeler;

import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.lattice.Folding;
import com.ariweiland.biophysics.lattice.Lattice;
import com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.peptide.Polypeptide;

import java.util.Queue;
import java.util.concurrent.PriorityBlockingQueue;

/**
 * @author Ari Weiland
 */
public class CurrentParallelModeler extends ParallelModeler {

    @Override
    protected PriorityBlockingQueue<Folding> initializeHeap(Polypeptide polypeptide) {
        PriorityBlockingQueue<Folding> initialHeap = new PriorityBlockingQueue<>();
        int size = polypeptide.size();
        // initialize the lattices
        Peptide first = polypeptide.get(0);
        Lattice line = new Lattice();
        line.put(0, 0, first);

        if (size > 1) {
            Peptide second = polypeptide.get(1);
            line.put(1, 0, second);

            // fill the queue initially.  this removes symmetrical solutions
            // if size == 2, the for loop will be ignored and none of this will matter
            double lowerBound = polypeptide.getMinEnergy()
                    - 2 * first.minInteraction()
                    + 3 * getFavorableWaterInteraction(first)
                    - 2 * second.minInteraction()
                    + 2 * getFavorableWaterInteraction(second);
            for (int i=2; i<size; i++) {
                Peptide next = polypeptide.get(i);
                Lattice bend = new Lattice(line);
                bend.put(i - 1, 1, next);
                line.put(i, 0, next);
                lowerBound += 2 * getFavorableWaterInteraction(next) - 2 * next.minInteraction();
                if (i == size - 1) {
                    lowerBound = bend.getEnergy();
                }
                initialHeap.add(new Folding(bend, i - 1, 1, i, lowerBound));
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
            double bound = folding.energyBound - 2 * p.minInteraction();
            for (Point.Direction d : Point.Direction.values()) {
                Point next = folding.lastPoint.getAdjacent(d);
                if (!folding.lattice.containsPoint(next)) {
                    Lattice l = new Lattice(folding.lattice);
                    l.put(next, p);
                    // though limiting the protein to the smallest possible rectangle is
                    // overly limiting, empirically it seems that limiting it to a rectangle
                    // of perimeter 4 larger does not seem to restrict the solution at all
                    if (l.boundingPerimeter() <= getPerimeterBound(polypeptide)) {
                        // subtract a water interaction where the next residue will end up
                        // note that if there is nowhere for the next residue, the foldings will be dropped on the next iteration
                        double nextBound = bound - getFavorableWaterInteraction(p);
                        if (nextIndex < size - 1) {
                            for (Point.Direction d1 : Point.Direction.values()) {
                                if (d1 != d.getReverse()) {
                                    if (l.containsPoint(next.getAdjacent(d1))) {
                                        Peptide adjacent = l.get(next.getAdjacent(d1));
                                        nextBound += p.interaction(adjacent) - getFavorableWaterInteraction(adjacent);
                                    } else {
                                        nextBound += getFavorableWaterInteraction(p);
                                    }
                                }
                            }
                        } else {
                            nextBound = l.getEnergy();
                        }
                        queue.add(new Folding(l, next, nextIndex, nextBound));
                    }
                }
            }
            return null;
        } else {
            return folding;
        }
    }
}
