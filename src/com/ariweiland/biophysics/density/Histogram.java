package com.ariweiland.biophysics.density;

import com.ariweiland.biophysics.peptide.Polypeptide;
import com.ariweiland.biophysics.peptide.Residue;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * @author Ari Weiland
 */
public abstract class Histogram {

    public abstract Map<Double, Integer> count(int dimension, Polypeptide polypeptide);

    public static void main(String[] args) {
        Histogram histogram = new PureMonteCarloHistogram(10000000);
//        Polypeptide polypeptide = new Polypeptide("(P)-(H)-(H)-(P)-(P)-(H)-(H)-(P)-(P)-(H)-(H)-(P)-(H)-(P)");
        Polypeptide polypeptide = new Polypeptide();
        for (int i=0; i<16; i++) {
            if (Math.random() < 0.4) {
                polypeptide.add(Residue.H);
            } else {
                polypeptide.add(Residue.P);
            }
        }
        System.out.println(polypeptide);
        System.out.println("Node count: " + polypeptide.size());
        System.out.println();
        long start = System.currentTimeMillis();
        Map<Double, Integer> count = histogram.count(3, polypeptide);
        long elapsed = System.currentTimeMillis() - start;
        System.out.println("Elapsed time: " + (elapsed / 1000.0) + " s");
        System.out.println("Bins : Counts");
        List<Double> keys = new ArrayList<>(count.keySet());
        Collections.sort(keys);
        for (double d : keys) {
            System.out.println(d + "\t" + count.get(d));
        }
    }

}
