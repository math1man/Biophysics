package com.ariweiland.biophysics.density;

import com.ariweiland.biophysics.peptide.Polypeptide;

import java.util.Map;

/**
 * @author Ari Weiland
 */
public interface Histogram {

    public Map<Double, Integer> count(int dimension, Polypeptide polypeptide);

}
