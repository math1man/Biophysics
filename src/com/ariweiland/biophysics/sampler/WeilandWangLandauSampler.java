package com.ariweiland.biophysics.sampler;

import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.RandomUtils;
import com.ariweiland.biophysics.lattice.MovableLattice;
import com.ariweiland.biophysics.lattice.PullMove;
import com.ariweiland.biophysics.lattice.RebridgeMove;
import com.ariweiland.biophysics.peptide.Polypeptide;

import java.util.List;
import java.util.Map;

/**
 * @author Ari Weiland
 */
public class WeilandWangLandauSampler extends WangLandauSampler {

    private int moveCount = 1;      // must be positive
    private double moveRatio = 0.2; // must be between 0 and 1 exclusive
    private boolean running;

    public WeilandWangLandauSampler() {
    }

    public WeilandWangLandauSampler(double flatness) {
        super(flatness);
    }

    public WeilandWangLandauSampler(double flatness, int moveCount) {
        super(flatness);
        this.moveCount = moveCount;
    }

    public WeilandWangLandauSampler(double flatness, int moveCount, double moveRatio) {
        super(flatness);
        this.moveCount = moveCount;
        this.moveRatio = moveRatio;
    }

    public int getMoveCount() {
        return moveCount;
    }

    public void setMoveCount(int moveCount) {
        this.moveCount = moveCount;
    }

    public double getMoveRatio() {
        return moveRatio;
    }

    public void setMoveRatio(double moveRatio) {
        this.moveRatio = moveRatio;
    }

    @Override
    public void terminate() {
        running = false;
    }

    @Override
    public Map<Double, Double> getDensity(int dimension, Polypeptide polypeptide) {
        running = true;
        int size = polypeptide.size();
        FValue f = new FValue(Math.E);
        g.clear();

        int count = 0;
        int pullCount = 0;
        int rebridgeCount = 0;
        while (Math.log(f.asDouble()) > F_FINAL && running) {
            h.clear();
            MovableLattice old = new MovableLattice(dimension, size);
            for (int i=0; i<size; i++) {
                old.put(new Point(i, 0, 0), polypeptide.get(i));
            }
            while (!isSufficientlyFlat() && running) {
                MovableLattice trial = new MovableLattice(old);
                List<PullMove> pullMoves = old.getPullMoves();
                int nOld = pullMoves.size();
                int pulls = 0;
                int rebridges = 0;
                for (int i=0; i<moveCount; i++) {
                    List<RebridgeMove> rebridgeMoves = trial.getRebridgeMoves();
                    if (rebridgeMoves.isEmpty() || RandomUtils.tryChance(moveRatio)) { // pull move
                        PullMove move = RandomUtils.selectRandom(pullMoves);
                        trial.pull(move);
                        pulls++;
                    } else {      // bond-rebridging move
                        RebridgeMove move = RandomUtils.selectRandom(rebridgeMoves);
                        trial.rebridge(move);
                        rebridges++;
                    }
                    pullMoves = trial.getPullMoves();
                }
                if (RandomUtils.tryChance(threshold(old, trial, ((double) nOld) / pullMoves.size()))) {
                    old = trial;
                    pullCount += pulls;
                    rebridgeCount += rebridges;
                }
                updateMaps(old.getEnergy(), f.asBigDecimal());
                count++;
                if (count % 1000000 == 0) {
                    System.out.println((count / 1000000) + "M trials");
                    System.out.println("\t" + (((double) pullCount)/count) + " pull move proportion");
                    System.out.println("\t" + (((double) rebridgeCount)/count) + " rebridge move proportion");
                }
            }
            reduceG();
            f.sqrt();
        }
        System.out.println(count + " total trials");
        System.out.println(pullCount + " total pull moves");
        System.out.println(rebridgeCount + " total rebridge moves");
        return convertG();
    }
}
