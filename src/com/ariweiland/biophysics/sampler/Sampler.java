package com.ariweiland.biophysics.sampler;

import com.ariweiland.biophysics.peptide.Polypeptide;

import java.util.*;

/**
 * @author Ari Weiland
 */
public abstract class Sampler {

    public abstract Map<Double, Double> getDensity(int dimension, Polypeptide polypeptide);
    
    public Map<Double, Double> normalize(Map<Double, Double> density) {
        double total = 0;
        for (double d : density.values()) {
            total += d;
        }
        Map<Double, Double> normalized = new HashMap<>();
        for (double d : density.keySet()) {
            normalized.put(d, density.get(d) / total);
        }
        return normalized;
    }

    public static void main(String[] args) {
        Sampler sampler = new WangLandauSampler(0.8, 10);
        Polypeptide polypeptide = new Polypeptide("(P)-(P)-(H)-(H)-(P)-(H)-(P)-(H)-(H)-(P)-(H)-(P)-(P)-(P)-(P)-(H)-(H)-(P)-(H)-(H)-(H)-(P)-(H)-(P)");
//        Polypeptide polypeptide = new Polypeptide();
//        for (int i=0; i<24; i++) {
//            if (RandomUtils.tryChance(0.4)) {
//                polypeptide.add(Residue.H);
//            } else {
//                polypeptide.add(Residue.P);
//            }
//        }
        System.out.println(polypeptide);
        System.out.println("Node count: " + polypeptide.size());
        System.out.println();
        long start = System.currentTimeMillis();
        Map<Double, Double> density = sampler.normalize(sampler.getDensity(2, polypeptide));
        long elapsed = System.currentTimeMillis() - start;
        System.out.println("Elapsed time: " + (elapsed / 1000.0) + " s");
        System.out.println("Bins\tCounts");
        List<Double> keys = new ArrayList<>(density.keySet());
        Collections.sort(keys);
        for (double d : keys) {
            System.out.println(d + " \t" + density.get(d));
        }
    }

}
