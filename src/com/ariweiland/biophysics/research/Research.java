package com.ariweiland.biophysics.research;

import com.ariweiland.biophysics.peptide.Polypeptide;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

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

    public static List<ErrorDataPoint> getErrorData(Map<Polypeptide, Map<String, Map<Double, Double>>> data,
                                                    ErrorFunction ef, String expectedFlag, String outcomeFlag) {
        List<ErrorDataPoint> dataPoints = new ArrayList<>();
        for (Polypeptide p : data.keySet()) {
            int size = p.size();
            Map<Double, Double> expected = data.get(p).get(expectedFlag);
            Map<Double, Double> outcome = data.get(p).get(outcomeFlag);

            if (outcome != null && outcome.size() == expected.size()) {
                dataPoints.add(new ErrorDataPoint(size, ef.error(expected, outcome)));
            }
        }
        Collections.sort(dataPoints);
        return dataPoints;
    }

    public static String asMathematicaCode(List<ErrorDataPoint> dataPoints) {
        StringBuilder dataString = new StringBuilder("{");
        for (ErrorDataPoint edp : dataPoints) {
            dataString.append("{").append(edp.length).append(",").append(edp.error).append("},\n");
        }
        // remove the final comma and newline, then add the closing }
        dataString.delete(dataString.length() - 2, dataString.length()).append("}");
        return dataString.toString();
    }

    public static void main(String[] args) {
        Map<Polypeptide, Map<String, Map<Double, Double>>> r30 = readFile("src/research/density/mass_r30.txt");
        Map<Polypeptide, Map<String, Map<Double, Double>>> r70 = readFile("src/research/density/mass_r70.txt");

        List<ErrorDataPoint> naiveDataPoints = getErrorData(r30, ErrorFunction.MEAN_SQUARED, "B", "N");
        List<ErrorDataPoint> wlDataPoints = getErrorData(r30, ErrorFunction.MEAN_SQUARED, "B", "W");

        System.out.println(asMathematicaCode(naiveDataPoints));
        System.out.println();
        System.out.println(asMathematicaCode(wlDataPoints));
    }
}
