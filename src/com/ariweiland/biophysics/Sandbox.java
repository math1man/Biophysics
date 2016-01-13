package com.ariweiland.biophysics.src.com.ariweiland.biophysics;

import com.ariweiland.biophysics.src.com.ariweiland.biophysics.lattice.SurfaceLattice;
import com.ariweiland.biophysics.src.com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.src.com.ariweiland.biophysics.peptide.Residue;

/**
 * @author Ari Weiland
 */
public class Sandbox {

    public static void main(String[] args) {
        SurfaceLattice lattice = new SurfaceLattice(Residue.H);
        lattice.put(0, 1, new Peptide(0, Residue.P));
        lattice.put(1, 1, new Peptide(1, Residue.H));
        lattice.put(2, 1, new Peptide(2, Residue.P));
        lattice.put(3, 1, new Peptide(3, Residue.H));
        lattice.put(3, 2, new Peptide(4, Residue.P));
        lattice.put(4, 2, new Peptide(5, Residue.P));
        lattice.put(4, 1, new Peptide(6, Residue.H));
        lattice.put(5, 1, new Peptide(7, Residue.P));
        lattice.put(6, 1, new Peptide(8, Residue.H));
        lattice.put(7, 1, new Peptide(9, Residue.P));
        lattice.put(8, 1, new Peptide(10, Residue.P));
        lattice.put(9, 1, new Peptide(11, Residue.P));
        lattice.visualize();
        System.out.println(lattice.getEnergy());
    }
}
