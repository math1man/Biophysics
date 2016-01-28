package com.ariweiland.biophysics;

import java.util.List;

/**
 * @author Ari Weiland
 */
public class RandomUtils {

    public static boolean tryChance(double fraction) {
        return fraction > 0 && (fraction >= 1 || Math.random() < fraction);
    }

    public static int randomInt(int max) {
        return (int) (Math.random() * max);
    }

    public static <T> T selectRandom(List<T> list) {
        if (list.isEmpty()) {
            throw new IllegalArgumentException("Nothing to select from");
        } else if (list.size() == 1) {
            return list.get(0);
        } else {
            return list.get(randomInt(list.size()));
        }
    }
}
