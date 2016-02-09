package com.ariweiland.biophysics.modeler;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.lattice.Folding;
import com.ariweiland.biophysics.lattice.BoundingLattice;
import com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.peptide.Polypeptide;

import java.util.Queue;
import java.util.concurrent.PriorityBlockingQueue;

/**
 * @author Ari Weiland
 */
public class CurrentParallelModeler extends ParallelModeler {

    public CurrentParallelModeler(int dimension) {
        super(dimension);
    }

    @Override
    protected PriorityBlockingQueue<Folding> initializeHeap(Polypeptide polypeptide) {
        PriorityBlockingQueue<Folding> initialHeap = new PriorityBlockingQueue<>(getSeedCount(polypeptide) * 4);
        int dim = getDimension();
        int size = polypeptide.size();
        // initialize the lattices
        Peptide first = polypeptide.get(0);
        BoundingLattice line = new BoundingLattice(dim, size);
        line.put(new Point(0, 0, 0), first);

        if (size > 1) {
            Peptide second = polypeptide.get(1);
            line.put(new Point(1, 0, 0), second);

            // fill the queue initially.  this removes symmetrical solutions
            // if size == 2, the for loop will be ignored and none of this will matter
            double lowerBound = polypeptide.getMinEnergy(dim)
                    - dim * first.minInteraction()
                    + (dim + 1) * getFavorableWaterInteraction(first)
                    - dim * second.minInteraction()
                    + dim * getFavorableWaterInteraction(second);
            for (int i=2; i<size; i++) {
                Peptide next = polypeptide.get(i);
                BoundingLattice bend = new BoundingLattice(line);
                Point point = new Point(i - 1, 1, 0);
                bend.put(point, next);
                line.put(new Point(i, 0, 0), next);
                lowerBound += (dim - 1) * 2 * getFavorableWaterInteraction(next) - (dim - 1) * 2 * next.minInteraction();
                if (i == size - 1) {
                    lowerBound = bend.getEnergy();
                }
                initialHeap.add(new Folding(bend, point, i, lowerBound));
            }
        }
        initialHeap.add(new Folding(line, new Point(size - 1, 0, 0), size - 1, line.getEnergy()));
        return initialHeap;
    }

    @Override
    public Folding iterate(Polypeptide polypeptide, Queue<Folding> queue) {
        int dim = getDimension();
        int size = polypeptide.size();
        Folding folding = queue.poll();
        int nextIndex = folding.index + 1;
        if (nextIndex < size) {
            Peptide p = polypeptide.get(nextIndex);
            for (Direction nextDir : Direction.values(dim)) {
                Point next = folding.lastPoint.getAdjacent(nextDir);
                if (!folding.lattice.containsPoint(next)) {
                    BoundingLattice l = new BoundingLattice(folding.lattice);
                    l.put(next, p);
                    // though limiting the protein to the smallest possible rectangle is
                    // overly limiting, empirically it seems that limiting it to a rectangle
                    // of perimeter 4 larger does not seem to restrict the solution at all
                    if (l.boundingPerimeter() <= getSurfaceBound(polypeptide)) {
                        // subtract a water interaction where the next residue will end up
                        // note that if there is nowhere for the next residue, the foldings will be dropped on the next iteration
                        double bound = folding.energyBound - (dim - 1) * 2 * p.minInteraction() - getFavorableWaterInteraction(p);
                        if (nextIndex < size - 1) {
                            for (Direction d : Direction.values(dim)) {
                                // the adjustments for the attached residue are already handled
                                if (d != nextDir.getReverse()) {
                                    if (l.containsPoint(next.getAdjacent(d))) {
                                        Peptide adjacent = l.get(next.getAdjacent(d));
                                        bound += p.interaction(adjacent) - getFavorableWaterInteraction(adjacent);
                                    } else {
                                        bound += getFavorableWaterInteraction(p);
                                    }
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
