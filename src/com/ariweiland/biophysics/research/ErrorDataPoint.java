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

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ErrorDataPoint that = (ErrorDataPoint) o;

        return Double.compare(that.error, error) == 0 && length == that.length;

    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = length;
        temp = Double.doubleToLongBits(error);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }
}
