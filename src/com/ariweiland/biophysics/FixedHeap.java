package com.ariweiland.biophysics;

import java.util.*;

/**
 * This is a heap implementation with a fixed size.
 * When the heap fills up to the fixed size, the next element replaces
 * the oldest leaf node. The ideal capacity of the heap is 2^n - 1 for
 * some integer n, and the implementation may not function properly if
 * the capacity is not of this form.
 *
 * @author Ari Weiland
 */
public class FixedHeap<T extends Comparable<T>> implements Queue<T> {

    private final Comparable<T>[] array;
    private int size;
    private int overflowAddIndex = 0;

    /**
     * Creates a hea with capacity 1023.
     */
    public FixedHeap() {
        this(1023);
    }

    /**
     * Creates a heap with the specified capacity.
     * Ideal capacity is of the form 2^n - 1.
     * @param capacity
     */
    public FixedHeap(int capacity) {
        this.array = new Comparable[capacity + 1]; // add in an extra index as a buffer
        this.size = 0;
    }

    @Override
    public boolean add(T t) {
        int currNode = size;
        int parentNode;
        boolean modified = false;
        if (size == array.length - 1) {
            // currNode refers to the buffer index at the end
            // parentNode refers to the next node to be replaced
            parentNode = overflowAddIndex + size / 2;
            overflowAddIndex = (overflowAddIndex + 1) % ((size + 1) / 2);
        } else {
            // currNode refers to the next available index
            // parentNode refers to the parent of that node
            parentNode = findParentNode(currNode);
            size++;
            modified = true;
        }
        array[currNode] = t;
        while (currNode > 0 && array[parentNode].compareTo(t) > 0) {
            array[currNode] = array[parentNode];
            array[parentNode] = t;
            currNode = parentNode;
            parentNode = findParentNode(currNode);
            modified = true;
        }
        return modified;
    }

    @Override
    public boolean remove(Object o) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean containsAll(Collection<?> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean addAll(Collection<? extends T> c) {
        boolean modified = false;
        for (T t : c) {
            modified = modified || add(t);
        }
        return modified;
    }

    @Override
    public boolean removeAll(Collection<?> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean retainAll(Collection<?> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean offer(T t) {
        return add(t);
    }

    @Override
    public T remove() {
        T t = poll();
        if (t == null) {
            throw new NoSuchElementException();
        } else {
            return t;
        }
    }

    @Override
    public T poll() {
        T largest = peek();
        size--;
        if (largest != null && size > 0) {
            T value = (T) array[size];
            array[0] = value;
            int currNode = 0;
            int bigChildNode = findBigChild(currNode);
            while (currNode < ((size) / 2) && array[bigChildNode].compareTo(value) < 0) {
                array[currNode] = array[bigChildNode];
                array[bigChildNode] = value;
                currNode = bigChildNode;
                bigChildNode = findBigChild(currNode);
            }
        }
        return largest;
    }

    @Override
    public T element() {
        if (size > 0) {
            return (T) array[0];
        } else {
            throw new NoSuchElementException();
        }
    }

    @Override
    public T peek() {
        if (size > 0) {
            return (T) array[0];
        } else {
            return null;
        }
    }

    @Override
    public int size() {
        return size;
    }

    @Override
    public boolean isEmpty() {
        return size == 0;
    }

    @Override
    public boolean contains(Object o) {
        return false;
    }

    @Override
    public void clear() {
        Arrays.fill(array, null);
        size = 0;
        overflowAddIndex = 0;
    }

    private int findParentNode(int index) {
        return (index + 1) / 2 - 1;
    }

    private int findLeftChildNode(int index) {
        return (index + 1) * 2 - 1;
    }

    private int findRightChildNode(int index) {
        return (index + 1) * 2;
    }

    private int findBigChild(int index) {
        int leftChildNode = findLeftChildNode(index);
        int rightChildNode = findRightChildNode(index);
        if (rightChildNode < size && array[leftChildNode].compareTo((T) array[rightChildNode]) > 0) {
            return rightChildNode;
        } else {
            return leftChildNode;
        }
    }

    private void buildHeap() {
        int half = size / 2 - 1;
        for (int r = half; r >= 0; r--) {
            int currNode = r;
            T value = (T) array[r];
            int bigChildNode = findBigChild(currNode);
            while (currNode <= half && array[bigChildNode].compareTo(value) < 0) {
                array[currNode] = array[bigChildNode];
                array[bigChildNode] = value;
                currNode = bigChildNode;
                bigChildNode = findBigChild(currNode);
            }
        }
    }

    /**
     * Returns the heap as a list in sorted order.
     * @return
     */
    public List<T> sort() {
        int oldSize = size;
        for (int i=size-1; i > 0; i--) {
            T next = poll();
            array[i] = next;
        }
        size = oldSize;
        List<T> list = new ArrayList<>(size);
        for (int i=0; i<size; i++) {
            list.add((T) array[i]);
        }
        buildHeap();
        return list;
    }

    @Override
    public Iterator<T> iterator() {
        return sort().iterator();
    }

    @Override
    public Object[] toArray() {
        return sort().toArray();
    }

    @Override
    public <T1> T1[] toArray(T1[] a) {
        return sort().toArray(a);
    }

    @Override
    public String toString() {
        Comparable<T>[] output = new Comparable[size];
        System.arraycopy(array, 0, output, 0, size);
        return Arrays.toString(output);
    }
}
