package com.ariweiland.biophysics.lattice;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.peptide.Residue;

/**
 * @author Ari Weiland
 */
public class FastLattice extends BoundingLattice {

    public FastLattice(int dimension) {
        super(dimension);
    }

    public FastLattice(int dimension, int initialCapacity) {
        super(dimension, initialCapacity);
    }

    public FastLattice(int dimension, Residue surface) {
        super(dimension, surface);
    }

    public FastLattice(int dimension, Residue surface, int initialCapacity) {
        super(dimension, surface, initialCapacity);
    }

    @Override
    public Peptide get(Point point) {
        if (hasSurface() && point.y == 0) {
            return new Peptide(-2, getSurface());
        } else {
            return lattice.get(point);
        }
    }

    @Override
    public void put(Point point, Peptide peptide) {
        for (Direction d : Direction.values(getDimension())) {
            Peptide adj = get(point.getAdjacent(d));
            if (adj != null) {
                // if they are not adjoining peptides
                if (adj.index != peptide.index + 1 && adj.index != peptide.index - 1) {
                    energy += peptide.interaction(adj);
                }
                energy -= adj.interaction(Residue.H2O);
            } else {
                energy += peptide.interaction(Residue.H2O);
            }
        }
        lattice.put(point, peptide);
    }

}
