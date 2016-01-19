package com.ariweiland.biophysics.density;

import com.ariweiland.biophysics.peptide.Polypeptide;

import java.util.Map;

/**
 * @author Ari Weiland
 */
public class BruteForceHistogram implements Histogram {

    @Override
    public Map<Double, Integer> count(int dimension, Polypeptide polypeptide) {
        int[] counter = new int[(dimension - 1) * 2 * polypeptide.size()];
        return null;
    }
}
