package com.ariweiland.biophysics.research;

import com.ariweiland.biophysics.peptide.Polypeptide;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Ari Weiland
 */
public class Research {

    public static Map<Polypeptide, Map<String, Map<Double, Double>>> readFile(String file) {
        BufferedReader br = null;
        List<String> lines = new ArrayList<>();
        try {
            br = new BufferedReader(new FileReader(file));
            String line = br.readLine();

            while (line != null) {
                lines.add(line);
                line = br.readLine();
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        } finally {
            if (br != null) {
                try {
                    br.close();
                } catch (IOException ignored) {}
            }
        }

        Map<Polypeptide, Map<String, Map<Double, Double>>> output = new HashMap<>();
        Polypeptide polypeptide = null;
        String type = null;
        for (String line : lines) {
            if (line.startsWith("Folding")) {
                polypeptide = new Polypeptide(line.substring(8));
                output.put(polypeptide, new HashMap<String, Map<Double, Double>>());
            } else if (line.startsWith("B") || line.startsWith("Na") || line.startsWith("W")) {
                type = line.substring(0, 1);
                output.get(polypeptide).put(type, new HashMap<Double, Double>());
            } else if (line.startsWith("-") || line.startsWith("0.0")) {
                String[] parts = line.split(" \t");
                output.get(polypeptide).get(type).put(Double.valueOf(parts[0]), Double.valueOf(parts[1]));
            }
        }
        return output;
    }

    public static double meanSquaredRatio(Map<Double, Double> expected, Map<Double, Double> outcome) {
        double sum = 0;
        for (double energy : expected.keySet()) {
            double e = expected.get(energy);
            double o = outcome.get(energy);
            double ratio = Math.log(o / e);
            sum += ratio * ratio;
        }
        return sum / expected.size();
    }

    public static double chiSquaredRatio(Map<Double, Double> expected, Map<Double, Double> outcome) {
        double sum = 0;
        for (double energy : expected.keySet()) {
            double e = expected.get(energy);
            double o = outcome.get(energy);
            double ratio = Math.log(o / e);
            sum += ratio * ratio / Math.log(e);
        }
        return -sum;
    }

    public static void main(String[] args) {
        Map<Polypeptide, Map<String, Map<Double, Double>>> r50 = readFile("src/research/density/mass_r50.txt");
        Map<Integer, List<Double>> naiveDataPoints = new HashMap<>();
        Map<Integer, List<Double>> wlDataPoints = new HashMap<>();
        for (Polypeptide p : r50.keySet()) {
            int size = p.size();
            if (!naiveDataPoints.containsKey(size)) {
                naiveDataPoints.put(size, new ArrayList<Double>());
                wlDataPoints.put(size, new ArrayList<Double>());
            }
            Map<Double, Double> expected = r50.get(p).get("B");
            Map<Double, Double> naive = r50.get(p).get("N");
            Map<Double, Double> wl = r50.get(p).get("W");

            if (naive != null && naive.size() == expected.size()) {
                naiveDataPoints.get(size).add(chiSquaredRatio(expected, naive));
            }
            if (wl != null && wl.size() == expected.size()) {
                wlDataPoints.get(size).add(chiSquaredRatio(expected, wl));
            }
        }
        
        StringBuilder naiveData = new StringBuilder("{");
        for (int size : naiveDataPoints.keySet()) {
            for (double error : naiveDataPoints.get(size)) {
                naiveData.append("{").append(size).append(",").append(error).append("},\n");
            }
        }
        naiveData.delete(naiveData.length() - 2, naiveData.length()).append("}");
        System.out.println(naiveData);
        System.out.println();

        StringBuilder wlData = new StringBuilder("{");
        for (int size : wlDataPoints.keySet()) {
            for (double error : wlDataPoints.get(size)) {
                wlData.append("{").append(size).append(",").append(error).append("},\n");
            }
        }
        wlData.delete(wlData.length() - 2, wlData.length()).append("}");
        System.out.println(wlData);
    }
}
