package com.ariweiland.biophysics.modeler;

import com.ariweiland.biophysics.lattice.Folding;
import com.ariweiland.biophysics.lattice.Lattice;
import com.ariweiland.biophysics.lattice.SurfaceLattice;
import com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.peptide.Polypeptide;
import com.ariweiland.biophysics.peptide.Residue;
import com.ariweiland.biophysics.Point;

import java.util.Queue;
import java.util.concurrent.PriorityBlockingQueue;

/**
 * @author Ari Weiland
 */
public class SurfaceModeler extends Modeler {

    private final Residue surface;

    public SurfaceModeler(int seedCount, Residue surface) {
        super(seedCount);
        this.surface = surface;
    }

    public Residue getSurface() {
        return surface;
    }

    public static int getMaxY(int n) {
        return getPerimBound(n) / 4 + 2;
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
                double lowerBound = polypeptide.getMinEnergy()
                        + getFavorableWaterInteraction(polypeptide.get(0))
                        - size * getAbsWaterSurfaceInteraction();
                for (k=0; k<=Math.abs(i-j) && k<size; k++) {
                    Peptide next = polypeptide.get(k);
                    int y;
                    if (i > j) {
                        y = i - k;
                    } else {
                        y = i + k;
                    }
                    lattice.put(0, y, next);
                    lowerBound += getFavorableWaterInteraction(next) - 2 * next.minInteraction() + getAbsWaterSurfaceInteraction();
                    if (y == 1) {
                        lowerBound += next.interaction(surface);
                    } else {
                        lowerBound += getFavorableWaterInteraction(next);
                    }
                }
                Point last;
                if (k < size) { // there is at least one residue left
                    Peptide next = polypeptide.get(k);
                    last = new Point(1, j);
                    lattice.put(last, next);
                    lowerBound += getFavorableWaterInteraction(next) - 2 * next.minInteraction() + getAbsWaterSurfaceInteraction();
                    if (j == 1) {
                        lowerBound += next.interaction(surface);
                    } else {
                        lowerBound += getFavorableWaterInteraction(next);
                    }
                } else {
                    last = new Point(0, j);
                }
                if (k > size - 2) { // the last residue added was the last residue
                    lowerBound = lattice.getEnergy();
                }
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

    private double getAbsWaterSurfaceInteraction() {
        return Math.abs(surface.interaction(Residue.H2O));
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
                if (!folding.lattice.containsPoint(next) && next.y < getMaxY(size)) {
                    SurfaceLattice l = new SurfaceLattice((SurfaceLattice) folding.lattice);
                    l.put(next, p);
                    // though limiting the protein to the smallest possible rectangle is
                    // overly limiting, empirically it seems that limiting it to a rectangle
                    // of perimeter 4 larger does not seem to restrict the solution at all
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
            return null;
        } else {
            return folding;
        }
    }
}
