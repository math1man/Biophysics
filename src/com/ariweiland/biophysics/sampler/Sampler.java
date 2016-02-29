package com.ariweiland.biophysics.sampler;

import com.ariweiland.biophysics.peptide.Polypeptide;
import com.ariweiland.biophysics.peptide.Residue;

import java.util.*;

/**
 * @author Ari Weiland
 */
public abstract class Sampler {

    public abstract void terminate();

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

    public static String asMathematicaCode(Map<Double, Double> density) {
        List<Double> keys = new ArrayList<>(density.keySet());
        Collections.sort(keys);
        StringBuilder dataString = new StringBuilder("{");
        for (double d : keys) {
            String cleaned = density.get(d).toString().replace("E", "*10^");
            dataString.append("{").append(d).append(",").append(cleaned).append("},\n");
        }
        // remove the final comma and newline, then add the closing }
        dataString.delete(dataString.length() - 2, dataString.length()).append("}");
        return dataString.toString();
    }

    public static void main(String[] args) {

        Sampler sampler = new BruteForceSurfaceSampler(Residue.H);
        Residue.setInteractionScheme(-7, -1, 0);
        Residue.setSurfaceInteractions(0, 0);
        int dimension = 2;

        Polypeptide polypeptide = Polypeptide.fibonacci(6);
//        Polypeptide polypeptide = Polypeptide.random(18, 0.5);
//        Polypeptide polypeptide = new Polypeptide("(H)-(P)-(P)-(P)-(H)-(P)-(H)-(P)-(H)-(H)-(P)-(P)-(P)-(H)-(P)-(H)-(P)-(H)");
        System.out.println(polypeptide);
        System.out.println("Node count: " + polypeptide.size());
        System.out.println();
        System.out.println(sampler.getClass().getSimpleName());
        long start = System.currentTimeMillis();
        Map<Double, Double> density = sampler.normalize(sampler.getDensity(dimension, polypeptide));
        long elapsed = System.currentTimeMillis() - start;
        System.out.println("Elapsed time: " + (elapsed / 1000.0) + " s");
        System.out.println("Bins\tCounts");
        List<Double> keys = new ArrayList<>(density.keySet());
        Collections.sort(keys);
        for (double d : keys) {
            System.out.println(d + " \t" + density.get(d));
        }
        System.out.println();
        System.out.println(asMathematicaCode(density));
    }
}
