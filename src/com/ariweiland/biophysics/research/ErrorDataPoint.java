package com.ariweiland.biophysics.research;

/**
 * @author Ari Weiland
 */
public class ErrorDataPoint implements Comparable<ErrorDataPoint> {

    public final int length;
    public final double error;

    public ErrorDataPoint(int length, double error) {
        this.length = length;
        this.error = error;
    }

    @Override
    public int compareTo(ErrorDataPoint o) {
        return this.length - o.length;
    }
}
