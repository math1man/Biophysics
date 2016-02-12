package com.ariweiland.biophysics.sampler;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.RandomUtils;
import com.ariweiland.biophysics.lattice.Lattice;
import com.ariweiland.biophysics.peptide.Polypeptide;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * This yields results pretty much identical to the basic naive Monte Carlo algorithm.
 * @author Ari Weiland
 */
public class NaiveWangLandauSampler extends WangLandauSampler {

    private boolean running;

    public NaiveWangLandauSampler() {}

    public NaiveWangLandauSampler(double flatness) {
        super(flatness);
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

        Lattice base = new Lattice(dimension, size);
        base.put(new Point(0, 0, 0), polypeptide.get(0));
        base.put(new Point(1, 0, 0), polypeptide.get(1));

        int count = 0;
        while (Math.log(f.asDouble()) > F_FINAL && running) {
            h.clear();
            Lattice old = null;
            while (!isSufficientlyFlat() && running) {
                Lattice trial = new Lattice(base);
                Point last = new Point(1, 0, 0);
                boolean isBoxedIn = false;
                for (int j=2; j<size && !isBoxedIn; j++) {
                    List<Direction> opens = new ArrayList<>();
                    for (Direction d : Direction.values(dimension)) {
                        if (!trial.contains(last.getAdjacent(d))) {
                            opens.add(d);
                        }
                    }
                    if (opens.isEmpty()) {
                        isBoxedIn = true;
                    } else {
                        Point next = last.getAdjacent(RandomUtils.selectRandom(opens));
                        trial.put(next, polypeptide.get(j));
                        last = next;
                    }
                }
                if (!isBoxedIn) {
                    if (old == null || RandomUtils.tryChance(threshold(old, trial, 1.0))) {
                        old = trial;
                    }
                    updateMaps(old.getEnergy(), f.asBigDecimal());
                    count++;
                    if (count % 1000000 == 0) {
                        System.out.println((count / 1000000) + "M trials");
                    }
                }
            }
            reduceG();
            f.sqrt();
        }
        System.out.println(count + " total states counted");
        return convertG();
    }
}
