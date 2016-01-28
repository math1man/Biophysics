package com.ariweiland.biophysics.lattice;

/**
 * @author Ari Weiland
 */
public class RebridgeMove {
    public final int i;
    public final int j;
    public final int k;
    public final int ip;
    public final int jp;
    public final int kp;

    /**
     * End or Parallel Rebridge Move
     * @param i
     * @param j
     * @param k
     */
    public RebridgeMove(int i, int j, int k) {
        this(i, j, k, -1, -1, -1);
    }

    /**
     * Antiparallel Rebridge Move
     * @param i
     * @param j
     * @param k
     * @param ip
     * @param jp
     * @param kp
     */
    public RebridgeMove(int i, int j, int k, int ip, int jp, int kp) {
        this.i = i;
        this.j = j;
        this.k = k;
        this.ip = ip;
        this.jp = jp;
        this.kp = kp;
    }

    @Override
    public String toString() {
        return "RebridgeMove{" +
                "i=" + i +
                ", j=" + j +
                ", k=" + k +
                ", ip=" + ip +
                ", jp=" + jp +
                ", kp=" + kp +
                '}';
    }
}
