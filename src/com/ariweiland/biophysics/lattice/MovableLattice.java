package com.ariweiland.biophysics.lattice;

import com.ariweiland.biophysics.Direction;
import com.ariweiland.biophysics.Point;
import com.ariweiland.biophysics.peptide.Peptide;
import com.ariweiland.biophysics.peptide.Residue;

import java.util.*;

/**
 * The standard lattice does not allow modifications to the lattice after a Residue has been dropped,
 * aside from the clear method. This lattice will have methods for both pull moves and bond-rebridging
 * moves. This lattice also does support surface size or bounding perimeter. It also dynamically
 * calculates energy instead of maintaining it throughout.
 * @author Ari Weiland
 */
public class MovableLattice extends Lattice {

    private final List<Point> pointSequence;
    private boolean modified = false;

    public MovableLattice(int dimension) {
        this(dimension, null);
    }

    public MovableLattice(int dimension, int initialCapacity) {
        this(dimension, null, initialCapacity);
    }

    public MovableLattice(int dimension, Residue surface) {
        super(dimension, surface);
        pointSequence = new ArrayList<>();
    }

    public MovableLattice(int dimension, Residue surface, int initialCapacity) {
        super(dimension, surface, initialCapacity);
        pointSequence = new ArrayList<>(initialCapacity);
    }

    // TODO this constructor is untested with normal Lattices
    public MovableLattice(Lattice lattice) {
        super(lattice);
        if (lattice instanceof MovableLattice) {
            pointSequence = new ArrayList<>(((MovableLattice) lattice).pointSequence);
        } else {
            pointSequence = new ArrayList<>(lattice.size());
            for (Point p : keySet()) {
                Peptide peptide = this.lattice.get(p);
                pointSequence.add(peptide.index, p);
            }
        }
    }

    @Override
    public void put(Point point, Peptide peptide) {
        if (getDimension() == 2 && point.z != 0) {
            throw new IllegalArgumentException("2D points cannot have a z-component");
        }
        if (hasSurface() && point.y < 1) {
            throw new IllegalArgumentException("Cannot put a point on or below the surface (y <= 0)");
        }
        if (containsPoint(point)) {
            throw new IllegalArgumentException("That point is already occupied");
        }
        for (Direction d : Direction.values(getDimension())) {
            Peptide adj = get(point.getAdjacent(d));
            if (adj != null) {
                // if they are not adjoining peptides
                if (adj.index != peptide.index + 1 && adj.index != peptide.index - 1) {
                    energy += peptide.interaction(adj);
                }
                energy -= adj.interaction(Residue.H2O);
            } else {
                energy += peptide.interaction(Residue.H2O);
            }
        }
        lattice.put(point, peptide);
        pointSequence.add(point);
    }

    @Override
    public void clear() {
        super.clear();
        pointSequence.clear();
    }

    @Override
    public double getEnergy() {
        if (modified) { // lazy update
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
        return energy;
    }

    @Override
    public int getSurfaceSize() {
        throw new UnsupportedOperationException();
    }

    @Override
    public int boundingPerimeter() {
        throw new UnsupportedOperationException();
    }

    /**
     * Applies a pull move on the peptide at index i to a random valid location.
     * Returns true if successful, or false if no valid location was found.
     * Throws an IllegalArgumentException if i specifies the end residue.
     * @param i
     * @return
     */
    public boolean pull(int i) {
        if (i == size() - 1) {
            throw new IllegalArgumentException("Cannot select the end residue");
        }
        Point point = pointSequence.get(i);
        Point next = pointSequence.get(i + 1);
        List<Direction> options = new ArrayList<>();
        for (Direction d : point.getDirectionTo(next).getNormals(getDimension())) {
            Point l = next.getAdjacent(d);
            Point c = point.getAdjacent(d);
            // position L is open and point c is open or occupied by peptide i-1
            if (!containsPoint(l) && (!containsPoint(c) || lattice.get(c).index == i - 1)) {
                options.add(d);
            }
        }
        if (options.isEmpty()) {
            return false;
        }
        Direction normal = selectRandom(options);
        Point l = next.getAdjacent(normal);
        Point c = point.getAdjacent(normal);
        Peptide peptide = lattice.remove(point);
        lattice.put(l, peptide); // first move point to L
        pointSequence.set(i, l);
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
        if (!containsPoint(c) && i > 0) {
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
        return true;
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
            List<Direction> options = new ArrayList<>();
            for (Direction d : Direction.values(getDimension())) {
                Point adj = point.getAdjacent(d);
                if (containsPoint(adj) && adj != pointSequence.get(i - 1)) {
                    options.add(d);
                }
            }
            if (options.isEmpty()) {
                return false;
            }
            Direction d = selectRandom(options);
            int j = lattice.get(point.getAdjacent(d)).index + 1;
            while (j < i) {
                Point p1 = pointSequence.get(j);
                Point p2 = pointSequence.get(i);
                Peptide temp = lattice.remove(p1);
                lattice.put(p1, lattice.remove(p2));
                lattice.put(p2, temp);
                pointSequence.set(i, p1);
                pointSequence.set(j, p2);
                j++;
                i--;
            }
            modified = true;
            return true;
        }
        Point next = pointSequence.get(i + 1);
        Direction direction = point.getDirectionTo(next);

        /*
        TODO:
        The following code in my opinion should be used to determine the normal direction, but the paper
        describes a different method, and initially I plan to follow the paper as closely as possible
        */

        List<Direction> options = new ArrayList<>();
        for (Direction d : point.getDirectionTo(next).getNormals(getDimension())) {
            if (!containsPoint(point.getAdjacent(d)) || !containsPoint(next.getAdjacent(d))) {
                return false;
            }
            int j = lattice.get(point.getAdjacent(d)).index;
            int k = lattice.get(next.getAdjacent(d)).index;
            // indices j and k are adjacent, and not consecutive to i and i+1
            if (Math.abs(j - k) == 1 && (k > i + 2 || j < i - 1)) {
                options.add(d);
            }
        }
        if (options.isEmpty()) {
            return false;
        }
        Direction normal = selectRandom(options);

//        Direction normal = selectRandom(Arrays.asList(direction.getNormals(getDimension())));

        if (!containsPoint(point.getAdjacent(normal)) || !containsPoint(next.getAdjacent(normal))) {
            return false;
        }

        int j = lattice.get(point.getAdjacent(normal)).index;
        int k = lattice.get(next.getAdjacent(normal)).index;
        if (k - j == 1) { // parallel case, type 2
            int m = j > i ? i + 1 : k;
            int n = j > i ? j : i;
            while (m < n) {
                Point p1 = pointSequence.get(m);
                Point p2 = pointSequence.get(n);
                Peptide temp = lattice.remove(p1);
                lattice.put(p1, lattice.remove(p2));
                lattice.put(p2, temp);
                pointSequence.set(n, p1);
                pointSequence.set(m, p2);
                m++;
                n--;
            }
        } else if (k - j == -1 && (k > i + 2 || j < i - 1)) { // antiparallel case, type 1
            // the second half of the above conditional prevents i, i+1, j, and k from being consecutive
            int loopStart;
            int loopEnd;
            if (i > j) { // loop is on the j side
                loopStart = j;
                loopEnd = i;
            } else {     // loop is on the k side
                loopStart = i + 1;
                loopEnd = k;
            }
            int ip = randomInt(loopEnd - loopStart) + loopStart;
            point = pointSequence.get(ip);
            next = pointSequence.get(ip + 1);
            direction = point.getDirectionTo(next);
            List<Direction> normals = Arrays.asList(direction.getNormals(getDimension()));
            Collections.shuffle(normals); // try all normals this time, but try them in a random order
            normal = null;
            int jp = 0;
            int kp = 0;
            for (int m=0; m<normals.size() && normal == null; m++) {
                Direction d = normals.get(m);
                if (!containsPoint(point.getAdjacent(d)) || !containsPoint(next.getAdjacent(d))) {
                    return false;
                }
                jp = lattice.get(point.getAdjacent(d)).index;
                kp = lattice.get(next.getAdjacent(d)).index;
                if (Math.abs(jp - kp) == 1 && ((kp > loopEnd && jp > loopEnd) || (kp < loopStart && jp < loopStart))) {
                    normal = d;
                }
            }
            if (normal == null) {
                return false;
            }
            int min = i;
            for (int a : Arrays.asList(i, k, ip, jp, kp)) {
                if (a < min) {
                    min = a;
                }
            }
            int max = i + 1;
            for (int a : Arrays.asList(i + 1, j, ip + 1, jp, kp)) {
                if (a > max) {
                    max = a;
                }
            }
            int[] changes = new int[max - min];
            // these booleans prevent the reordering to try to go both ways on each swap
            boolean ij = true;
            boolean ik = true;
            boolean ijp = true;
            boolean ikp = true;
            int change = 1;
            // changes[n] keeps track of the point/peptide reordering
            // the (min+n)th peptide should be moved to the changes[n]'th point of the original sequence
            changes[0] = min;
            for (int n=1; n<changes.length; n++) {
                int index = n + min;            // index specifies which peptide is being moved
                int lastPoint = changes[n - 1]; // lastPoint specifies the previous point in the modified sequence
                int nextPoint;                  // nextPoint will specify the point to follow lastPoint
                if (lastPoint == i && ij) {
                    nextPoint = j;
                    ij = false;
                    change = 1;
                } else if (lastPoint == j && ij) {
                    nextPoint = i;
                    ij = false;
                    change = -1;
                } else if (lastPoint == i + 1 && ik) {
                    nextPoint = k;
                    ik = false;
                    change = -1;
                } else if (lastPoint == k && ik) {
                    nextPoint = i + 1;
                    ik = false;
                    change = 1;
                } else if (lastPoint == ip && ijp) {
                    nextPoint = jp;
                    ijp = false;
                    change = jp - kp;
                } else if (lastPoint == jp && ijp) {
                    nextPoint = ip;
                    ijp = false;
                    change = -1;
                } else if (lastPoint == ip + 1 && ikp) {
                    nextPoint = kp;
                    ikp = false;
                    change = kp - jp;
                } else if (lastPoint == kp && ikp) {
                    nextPoint = ip + 1;
                    ikp = false;
                    change = 1;
                } else {
                    nextPoint = lastPoint + change;
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
        } else { // not a valid selection
            return false;
        }
        modified = true;
        return true;
    }

    private static int randomInt(int max) {
        return (int) (Math.random() * max);
    }

    private static <T> T selectRandom(List<T> list) {
        if (list.isEmpty()) {
            throw new IllegalArgumentException("Nothing to select from");
        } else if (list.size() == 1) {
            return list.get(0);
        } else {
            return list.get(randomInt(list.size()));
        }
    }

}
