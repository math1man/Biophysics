package com.ariweiland.biophysics.modeler;

import com.ariweiland.biophysics.*;
import com.ariweiland.biophysics.lattice.Folding;
import com.ariweiland.biophysics.lattice.Lattice;
import com.ariweiland.biophysics.lattice.SurfaceLattice;
import com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.peptide.Polypeptide;
import com.ariweiland.biophysics.peptide.Residue;

import java.util.Queue;
import java.util.concurrent.PriorityBlockingQueue;

/**
 * @author Ari Weiland
 */
public class SurfaceModeler extends ParallelModeler {

    private final Residue surface;

    public SurfaceModeler(int seedCount, Residue surface) {
        super(seedCount);
        this.surface = surface;
    }

    public Residue getSurface() {
        return surface;
    }

    public static int getMaxY(int n) {
        return getPerimBound(n) / 2 + 1;
    }

    public Lattice fold(Polypeptide polypeptide) {
        // fill the queue initially.  this avoids symmetrical solutions
        PriorityBlockingQueue<Folding> initialHeap = new PriorityBlockingQueue<>();
        int size = polypeptide.size();
        int maxY = getMaxY(size);
        for (int i=1; i<maxY; i++) {
            for (int j=1; j<maxY; j++) {
                SurfaceLattice lattice = new SurfaceLattice(surface);
                int k;
                double lowerBound = polypeptide.getMinEnergy();
                for (k=0; k<=Math.abs(i-j); k++) {
                    Peptide next = polypeptide.get(k);
                    if (i > j) {
                        lattice.put(0, i - k, next);
                    } else {
                        lattice.put(0, i + k, next);
                    }
                    lowerBound -= 2 * next.minInteraction();
                }
                Peptide next = polypeptide.get(k);
                Point last = new Point(1, j);
                lattice.put(last, next);
                lowerBound -= 2 * next.minInteraction();
                initialHeap.add(new Folding(lattice, last, k, lowerBound));
            }
        }

        // iterate a few times to make the initial heap bigger
        int count = 0;
        while (count < getSeedCount()) {
            Folding solution = iterate(polypeptide, initialHeap);
            if (solution != null) {
                return solution.lattice;
            }
            count++;
        }

        int processors = Runtime.getRuntime().availableProcessors();
        return parallelize(polypeptide, initialHeap, processors, MAX_HEAP_SIZE / processors - 1);
    }

    // TODO: this should be good to go!
    public Folding iterate(Polypeptide polypeptide, Queue<Folding> queue) {
        int size = polypeptide.size();
        Folding folding = queue.poll();
        int nextIndex = folding.index + 1;
        if (nextIndex < size) {
            Peptide p = polypeptide.get(nextIndex);
            double bound = folding.energyBound - 2 * p.minInteraction();
            for (Point.Direction d : Point.Direction.values()) {
                Point next = folding.lastPoint.getAdjacent(d);
                if (!folding.lattice.containsPoint(next) && next.y > 0 && next.y < getMaxY(size)) {
                    Lattice l = new Lattice(folding.lattice);
                    l.put(next, p);
                    // though limiting the protein to the smallest possible rectangle is
                    // overly limiting, empirically it seems that limiting it to a rectangle
                    // of perimeter 4 larger does not seem to restrict the solution at all
                    double lb;
                    if (nextIndex < size - 1) {
                        lb = bound - l.get(next.getAdjacent(d.getReverse())).interaction(Residue.H2O);
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
                    queue.add(new Folding(l, next, nextIndex, lb));
                }
            }
            return null;
        } else {
            return folding;
        }
    }
}
