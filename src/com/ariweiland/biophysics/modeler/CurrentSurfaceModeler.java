package com.ariweiland.biophysics.modeler;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.lattice.Folding;
import com.ariweiland.biophysics.lattice.SurfaceLattice;
import com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.peptide.Polypeptide;
import com.ariweiland.biophysics.peptide.Residue;

import java.util.Queue;

/**
 * @author Ari Weiland
 */
public class CurrentSurfaceModeler extends SurfaceModeler {


    public CurrentSurfaceModeler(Residue surface) {
        super(surface);
    }

    /**
     * This helper method returns the min interaction of the surface,
     * adjusted such that the water-surface interaction is treated as zero.
     * @return
     */
    private double getAdjustedSurfaceMinInteraction() {
        return getSurface().minInteraction() - getSurface().interaction(Residue.H2O);
    }

    protected double getInitialEnergyBound(Polypeptide polypeptide) {
        return polypeptide.getMinEnergy()
                + getFavorableWaterInteraction(polypeptide.get(0))
                + polypeptide.size() * getAdjustedSurfaceMinInteraction();
    }

    /**
     * This is a helper method that calculates the modification to the lower energy bound
     * in the initial seeding stage of filling the heap. It therefore assumes that no
     * peptide is adjacent to any other peptide in the polypeptide, though they may be
     * adjacent to the surface.
     *
     * First, it adds at least one of the peptide's favorable water interaction, and
     * subtracts the peptide's min interactions. Then, if it is adjacent to the surface,
     * it adds the interaction with the surface and subtracts any favorable water-surface
     * interaction. Otherwise, it adds another of the peptide's favorable water interactions,
     * and adds any unfavorable water-surface interaction because this peptide cannot block
     * one of those.
     *
     * @param y
     * @param p
     * @return
     */
    @Override
    protected double getBoundAdjust(int y, Peptide p) {
        double boundAdjust = getFavorableWaterInteraction(p) - 2 * p.minInteraction();
        if (y == 1) {
            boundAdjust += p.interaction(getSurface()) - getAdjustedSurfaceMinInteraction();
        } else {
            boundAdjust += getFavorableWaterInteraction(p);
        }
        return boundAdjust;
    }

    @Override
    public Folding iterate(Polypeptide polypeptide, Queue<Folding> queue) {
        int size = polypeptide.size();
        Folding folding = queue.poll();
        int nextIndex = folding.index + 1;
        if (nextIndex < size) {
            Peptide p = polypeptide.get(nextIndex);
            // try to add the peptide in every direction
            for (Direction d : Direction.values()) {
                Point next = folding.lastPoint.getAdjacent(d);
                if (!folding.lattice.containsPoint(next) && next.y < getMaxY(polypeptide)) {
                    SurfaceLattice l = new SurfaceLattice((SurfaceLattice) folding.lattice);
                    l.put(next, p);
                    // set the bound from the previous bound, minus the min interactions for this peptide,
                    // minus one favorable water interaction which
                    double bound = folding.energyBound - 2 * p.minInteraction() - getFavorableWaterInteraction(p);
                    if (nextIndex < size - 1) {
                        for (Direction d1 : Direction.values()) {
                            if (d1 != d.getReverse()) {
                                if (l.containsPoint(next.getAdjacent(d1))) {
                                    Peptide adjacent = l.get(next.getAdjacent(d1));
                                    bound += p.interaction(adjacent) - getFavorableWaterInteraction(adjacent);
                                } else {
                                    bound += getFavorableWaterInteraction(p);
                                }
                            }
                        }
                        if (next.y > 1) {
                            bound -= getAdjustedSurfaceMinInteraction();
                        }
                    } else {
                        bound = l.getEnergy();
                    }
                    queue.add(new Folding(l, next, nextIndex, bound));
                }
            }
            return null;
        } else {
            return folding;
        }
    }
}
