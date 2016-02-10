package com.ariweiland.biophysics.lattice;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.RandomUtils;
import com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.peptide.Residue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * The standard lattice does not allow modifications to the lattice after a Residue has been dropped,
 * aside from the clear method. This lattice will have methods for both pull moves and bond-rebridging
 * moves. This lattice also does support surface size or bounding perimeter. It also dynamically
 * calculates energy instead of maintaining it throughout.
 * @author Ari Weiland
 */
public class MovableLattice extends Lattice {

    private final List<Point> pointSequence;

    public MovableLattice(int dimension) {
        this(dimension, null);
    }

    public MovableLattice(int dimension, Residue surface) {
        super(dimension, surface);
        pointSequence = new ArrayList<>();
    }

    public MovableLattice(int dimension, int initialCapacity) {
        this(dimension, initialCapacity, null);
    }

    public MovableLattice(int dimension, int initialCapacity, Residue surface) {
        super(dimension, initialCapacity, surface);
        pointSequence = new ArrayList<>(initialCapacity);
    }

    public MovableLattice(MovableLattice lattice) {
        super(lattice);
        this.pointSequence = new ArrayList<>(lattice.pointSequence);
    }

    @Override
    public void put(Point point, Peptide peptide) {
        super.put(point, peptide);
        pointSequence.add(point);
    }

    @Override
    public void clear() {
        super.clear();
        pointSequence.clear();
    }

    /**
     * TODO: make pull moves surface-safe
     * @return
     */
    public List<PullMove> getPullMoves() {
        List<PullMove> moves = new ArrayList<>();
        for (int i=0; i<size()-1; i++) {
            Point point = pointSequence.get(i);
            Point next = pointSequence.get(i + 1);
            for (Direction d : point.getDirectionTo(next).getNormals(getDimension())) {
                Point l = next.getAdjacent(d);
                Point c = point.getAdjacent(d);
                // position L is open and i is 0, point c is open, or it's occupied by peptide i-1
                if (!contains(l) && (i == 0 || !lattice.containsKey(c) || lattice.get(c).index == i - 1)) {
                    moves.add(new PullMove(i, d)); // TODO try reversing the above &&
                }
            }
        }
        return moves;
    }

    /**
     * Applies a pull move as specified by the input parameter.
     * @param move
     */
    public void pull(PullMove move) {
        int i = move.index;
        Point point = pointSequence.get(i);
        Point next = pointSequence.get(i + 1);
        Direction normal = move.direction;
        Point l = next.getAdjacent(normal);
        Point c = point.getAdjacent(normal);
        Peptide peptide = lattice.remove(point);
        lattice.put(l, peptide); // first move point to L
        pointSequence.set(i, l);
        // update energy
        for (Direction d : Direction.values(getDimension())) {
            Point adj1 = point.getAdjacent(d);
            if (!adj1.equals(next) && !adj1.equals(c)) {
                energy -= peptide.interaction(lattice.get(adj1));
            }
            Point adj2 = l.getAdjacent(d);
            if (!adj2.equals(next) && !adj2.equals(c)) {
                energy += peptide.interaction(lattice.get(adj2));
            }
        }
        // pull the move along
        if (i > 0 && !lattice.containsKey(c)) {
            int j = i - 1;
            next = l;           // the new location of point j + 1 (used to check adjacency)
            Point two = c;      // the former location of point j + 2 and the new location of point j
            Point one = point;  // the former location of point j + 1 (temporary variable)
            while (j >= 0) {
                point = pointSequence.get(j); // get the next point
                if (!point.isAdjacentTo(next)) {
                    peptide = lattice.remove(point);
                    lattice.put(two, peptide); // move point to two
                    pointSequence.set(j, two);
                    // update energy
                    for (Direction d : Direction.values(getDimension())) {
                        Point adj1 = point.getAdjacent(d);
                        if (!adj1.equals(one)) {
                            energy -= peptide.interaction(lattice.get(adj1));
                        }
                        Point adj2 = two.getAdjacent(d);
                        if (!adj2.equals(one)) {
                            energy += peptide.interaction(lattice.get(adj2));
                        }
                    }
                    next = two;
                    two = one;
                    one = point;
                    j--;
                } else {
                    j = -1; // break
                }
            }
        }
    }

    public List<RebridgeMove> getRebridgeMoves() {
        List<RebridgeMove> moves = new ArrayList<>();
        int last = size() - 1;
        for (int i=0; i<last; i++) {
            Point point = pointSequence.get(i);
            Point next = pointSequence.get(i + 1);
            for (Direction d : point.getDirectionTo(next).getNormals(getDimension())) {
                Point pj = point.getAdjacent(d);
                Point pk = next.getAdjacent(d);
                if (lattice.containsKey(pj) && lattice.containsKey(pk)) {
                    int j = lattice.get(pj).index;
                    int k = lattice.get(pk).index;
                    // indices j and k are adjacent, and not consecutive to i and i+1
                    if (Math.abs(j - k) == 1) {
                        if (k > j) { // parallel
                            moves.add(new RebridgeMove(i, j, k));
                        } else if (k > i + 2 || j < i - 1) { // antiparallel
                            // need to make sure that i, i+1, j, and k are not consecutive
                            int loopStart;
                            int loopEnd;
                            if (i > j) { // loop is on the j side
                                loopStart = j;
                                loopEnd = i;
                            } else {     // loop is on the k side
                                loopStart = i + 1;
                                loopEnd = k;
                            }
                            for (int ip = loopStart; ip < loopEnd; ip++) {
                                Point pointP = pointSequence.get(ip);
                                Point nextP = pointSequence.get(ip + 1);
                                for (Direction d2 : pointP.getDirectionTo(nextP).getNormals(getDimension())) {
                                    Point pjp = pointP.getAdjacent(d2);
                                    Point pkp = nextP.getAdjacent(d2);
                                    if (lattice.containsKey(pjp) && lattice.containsKey(pkp)) {
                                        int jp = lattice.get(pjp).index;
                                        int kp = lattice.get(pkp).index;
                                        if (Math.abs(jp - kp) == 1 && ((kp > loopEnd && jp > loopEnd) || (kp < loopStart && jp < loopStart))) {
                                            moves.add(new RebridgeMove(i, j, k, ip, jp, kp));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        // end rebridge moves
//        Point point = pointSequence.get(last);
//        for (Direction d : Direction.values(getDimension())) {
//            Point adj = point.getAdjacent(d);
//            if (!pointSequence.get(last - 1).equals(adj) && lattice.containsKey(adj)) {
//                int j = lattice.get(point.getAdjacent(d)).index;
//                moves.add(new RebridgeMove(last, j, j + 1));
//            }
//        }
        return moves;
    }

    /**
     * Applies a rebridge move as specified by the input parameter.
     * @param move
     */
    public void rebridge(RebridgeMove move) {
        int i = move.i;
        int j = move.j;
        int k = move.k;
        int ip = move.ip;
        if (ip == -1) { // end case or parallel case, type 2
            reversePeptideSequence(j > i ? i + 1 : k, j > i ? j : i);
        } else { // antiparallel case, type 1
            int jp = move.jp;
            int kp = move.kp;
            int min = i;
            for (int a : Arrays.asList(k, ip, jp, kp)) {
                if (a < min) {
                    min = a;
                }
            }
            int max = i + 1;
            for (int a : Arrays.asList(j, ip + 1, jp, kp)) {
                if (a > max) {
                    max = a;
                }
            }
            int[] changes = new int[max - min];
            // these booleans prevent the reordering to try to go both ways on each swap
            boolean ij = true; boolean ik = true; boolean ijp = true; boolean ikp = true;
            int direction = 1; // direction indicates whether we should increment or decrement
            // changes[n] keeps track of the point/peptide reordering
            // the (min+n)th peptide should be moved to the changes[n]'th point of the original sequence
            changes[0] = min;
            for (int n = 1; n < changes.length; n++) {
                int index = n + min;            // index specifies which peptide is being moved
                int lastPoint = changes[n - 1]; // lastPoint specifies the previous point in the modified sequence
                int nextPoint;                  // nextPoint will specify the point to follow lastPoint
                if (lastPoint == i && ij) {
                    nextPoint = j;
                    ij = false;
                    direction = 1;
                } else if (lastPoint == j && ij) {
                    nextPoint = i;
                    ij = false;
                    direction = -1;
                } else if (lastPoint == i + 1 && ik) {
                    nextPoint = k;
                    ik = false;
                    direction = -1;
                } else if (lastPoint == k && ik) {
                    nextPoint = i + 1;
                    ik = false;
                    direction = 1;
                } else if (lastPoint == ip && ijp) {
                    nextPoint = jp;
                    ijp = false;
                    direction = jp - kp;
                } else if (lastPoint == jp && ijp) {
                    nextPoint = ip;
                    ijp = false;
                    direction = -1;
                } else if (lastPoint == ip + 1 && ikp) {
                    nextPoint = kp;
                    ikp = false;
                    direction = kp - jp;
                } else if (lastPoint == kp && ikp) {
                    nextPoint = ip + 1;
                    ikp = false;
                    direction = 1;
                } else {
                    nextPoint = lastPoint + direction;
                }
                changes[n] = nextPoint;
                Point oldPoint = pointSequence.get(index);    // oldPoint is the location of the index'th peptide
                int newIndex = nextPoint;                     // newIndex is the current location of the point
                while (newIndex < index) {                    // where the index'th peptide needs to go
                    newIndex = changes[newIndex - min];       // this while loop handles complex swapping behavior
                }
                Point newPoint = pointSequence.get(newIndex); // newPoint is the point that peptide needs to be moved to
                if (!oldPoint.equals(newPoint)) {             // if newPoint and oldPoint are different, do the swap
                    Peptide temp = lattice.remove(oldPoint);
                    lattice.put(oldPoint, lattice.remove(newPoint));
                    lattice.put(newPoint, temp);
                    pointSequence.set(newIndex, oldPoint);
                    pointSequence.set(index, newPoint);
                }
            }
        }
        energy = 0;
        for (Point p : pointSequence) {
            Peptide peptide = lattice.get(p);
            for (Direction d : Direction.values(getDimension())) {
                Peptide adj = lattice.get(p.getAdjacent(d));
                // don't count if lower index to avoid double counting
                // don't count if the next index because that is not an interaction
                if (adj == null || adj.index > peptide.index + 1) {
                    energy += peptide.interaction(adj);
                }
            }
        }
    }

    /**
     * Applies a bond-rebridging move on the peptide at index i to a random valid location.
     * Returns true if successful, or false if no valid location was found.
     * Throws an IllegalArgumentException if i specifies the end residue.
     * @param i
     * @return
     */
    public boolean rebridge(int i) {
        Point point = pointSequence.get(i);
        if (i == size() - 1) { // the unique end rebridging case
            // pick a direction to rebridge
            List<Direction> options = new ArrayList<>();
            for (Direction d : Direction.values(getDimension())) {
                Point adj = point.getAdjacent(d);
                if (contains(adj) && adj != pointSequence.get(i - 1)) {
                    options.add(d);
                }
            }
            // if no directions, stop
            if (options.isEmpty()) {
                return false;
            }
            // fix the sequence by reversing the order of peptides between indexes j+1 and i
            reversePeptideSequence(lattice.get(point.getAdjacent(RandomUtils.selectRandom(options))).index + 1, i);
        } else {
            Point next = pointSequence.get(i + 1);
            // pick a direction to rebridge
            List<Direction> options = new ArrayList<>();
            for (Direction d : point.getDirectionTo(next).getNormals(getDimension())) {
                Point pj = point.getAdjacent(d);
                Point pk = next.getAdjacent(d);
                if (lattice.containsKey(pj) && lattice.containsKey(pk)) {
                    int j = lattice.get(pj).index;
                    int k = lattice.get(pk).index;
                    // indices j and k are adjacent, and not consecutive to i and i+1
                    if (Math.abs(j - k) == 1 && (k > j || k > i + 2 || j < i - 1)) {
                        options.add(d);
                    }
                }
            }
            // if no directions, stop
            if (options.isEmpty()) {
                return false;
            }
            // begin rebridging
            Direction normal = RandomUtils.selectRandom(options);
            int j = lattice.get(point.getAdjacent(normal)).index;
            int k = lattice.get(next.getAdjacent(normal)).index;
            if (k - j == 1) { // parallel case, type 2
                // in this case, we need only reverse the sequence between the rebridged bonds
                if (i > j) {
                    reversePeptideSequence(k, i);
                } else {
                    reversePeptideSequence(i + 1, j);
                }
            } else { // antiparallel case, type 1
                int loopStart;
                int loopEnd;
                if (i > j) { // loop is on the j side
                    loopStart = j;
                    loopEnd = i;
                } else {     // loop is on the k side
                    loopStart = i + 1;
                    loopEnd = k;
                }
                int ip = RandomUtils.randomInt(loopEnd - loopStart) + loopStart;
                Point pointP = pointSequence.get(ip);
                Point nextP = pointSequence.get(ip + 1);
                List<Direction> optionsP = new ArrayList<>();
                for (Direction d : point.getDirectionTo(next).getNormals(getDimension())) {
                    Point pjp = pointP.getAdjacent(d);
                    Point pkp = nextP.getAdjacent(d);
                    if (lattice.containsKey(pjp) && lattice.containsKey(pkp)) {
                        int jp = lattice.get(pjp).index;
                        int kp = lattice.get(pkp).index;
                        if (Math.abs(jp - kp) == 1 && ((kp > loopEnd && jp > loopEnd) || (kp < loopStart && jp < loopStart))) {
                            optionsP.add(d);
                        }
                    }
                }
                if (optionsP.isEmpty()) {
                    return false;
                }
                Direction normalP = RandomUtils.selectRandom(optionsP);
                int jp = lattice.get(pointP.getAdjacent(normalP)).index;
                int kp = lattice.get(nextP.getAdjacent(normalP)).index;
                int min = i;
                for (int a : Arrays.asList(k, ip, jp, kp)) {
                    if (a < min) {
                        min = a;
                    }
                }
                int max = i + 1;
                for (int a : Arrays.asList(j, ip + 1, jp, kp)) {
                    if (a > max) {
                        max = a;
                    }
                }
                int[] changes = new int[max - min];
                // these booleans prevent the reordering to try to go both ways on each swap
                boolean ij = true; boolean ik = true; boolean ijp = true; boolean ikp = true;
                int direction = 1; // direction indicates whether we should increment or decrement
                // changes[n] keeps track of the point/peptide reordering
                // the (min+n)th peptide should be moved to the changes[n]'th point of the original sequence
                changes[0] = min;
                for (int n = 1; n < changes.length; n++) {
                    int index = n + min;            // index specifies which peptide is being moved
                    int lastPoint = changes[n - 1]; // lastPoint specifies the previous point in the modified sequence
                    int nextPoint;                  // nextPoint will specify the point to follow lastPoint
                    if (lastPoint == i && ij) {
                        nextPoint = j;
                        ij = false;
                        direction = 1;
                    } else if (lastPoint == j && ij) {
                        nextPoint = i;
                        ij = false;
                        direction = -1;
                    } else if (lastPoint == i + 1 && ik) {
                        nextPoint = k;
                        ik = false;
                        direction = -1;
                    } else if (lastPoint == k && ik) {
                        nextPoint = i + 1;
                        ik = false;
                        direction = 1;
                    } else if (lastPoint == ip && ijp) {
                        nextPoint = jp;
                        ijp = false;
                        direction = jp - kp;
                    } else if (lastPoint == jp && ijp) {
                        nextPoint = ip;
                        ijp = false;
                        direction = -1;
                    } else if (lastPoint == ip + 1 && ikp) {
                        nextPoint = kp;
                        ikp = false;
                        direction = kp - jp;
                    } else if (lastPoint == kp && ikp) {
                        nextPoint = ip + 1;
                        ikp = false;
                        direction = 1;
                    } else {
                        nextPoint = lastPoint + direction;
                    }
                    changes[n] = nextPoint;
                    Point oldPoint = pointSequence.get(index);    // oldPoint is the location of the index'th peptide
                    int newIndex = nextPoint;                     // newIndex is the current location of the point
                    while (newIndex < index) {                    // where the index'th peptide needs to go
                        newIndex = changes[newIndex - min];       // this while loop handles complex swapping behavior
                    }
                    Point newPoint = pointSequence.get(newIndex); // newPoint is the point that peptide needs to be moved to
                    if (!oldPoint.equals(newPoint)) {             // if newPoint and oldPoint are different, do the swap
                        Peptide temp = lattice.remove(oldPoint);
                        lattice.put(oldPoint, lattice.remove(newPoint));
                        lattice.put(newPoint, temp);
                        pointSequence.set(newIndex, oldPoint);
                        pointSequence.set(index, newPoint);
                    }
                }
            }
        }
        // update energy
        energy = 0;
        for (Point p : pointSequence) {
            Peptide peptide = lattice.get(p);
            for (Direction d : Direction.values(getDimension())) {
                Peptide adj = lattice.get(p.getAdjacent(d));
                // don't count if lower index to avoid double counting
                // don't count if the next index because that is not an interaction
                if (adj == null || adj.index > peptide.index + 1) {
                    energy += peptide.interaction(adj);
                }
            }
        }
        return true;
    }

    private void reversePeptideSequence(int firstIndex, int lastIndex) {
        while (firstIndex < lastIndex) {
            // get the points
            Point p1 = pointSequence.get(firstIndex);
            Point p2 = pointSequence.get(lastIndex);
            // swap the peptides
            Peptide temp = lattice.remove(p1);
            lattice.put(p1, lattice.remove(p2));
            lattice.put(p2, temp);
            // swap the points in the sequence
            pointSequence.set(lastIndex, p1);
            pointSequence.set(firstIndex, p2);
            // move to next
            firstIndex++;
            lastIndex--;
        }
    }

}
