package com.ariweiland.biophysics;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/**
 * This is a heap implementation with a fixed size.
 * When the heap fills up to the fixed size, the next element
 * @author Ari Weiland
 */
public class FixedHeap<T extends Comparable<T>> implements Iterable<T> {

    private Comparable<T>[] array;
    private int size;
    private int overflowAddIndex = 0;

    public FixedHeap() {
        this(1023);
    }

    public FixedHeap(int initialCap) {
        array = new Comparable[initialCap + 1]; // add in an extra index as a buffer
        size = 0;
    }

    public FixedHeap(T... array) {
        this.array = new Comparable[array.length + 1]; // add in an extra place as a buffer
        System.arraycopy(array, 0, this.array, 0, array.length);
        size = array.length;
        buildHeap();
    }

    public void push(T t) {
        int currNode = size;
        int parentNode;
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
        }
        array[currNode] = t;
        while (currNode > 0 && array[parentNode].compareTo(t) > 0) {
            array[currNode] = array[parentNode];
            array[parentNode] = t;
            currNode = parentNode;
            parentNode = findParentNode(currNode);
        }
    }

    public T peek() {
        if (size > 0) {
            return (T) array[0];
        } else {
            return null;
        }
    }

    public T pop() {
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

    public int size() {
        return size;
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

    public List<T> sort() {
        int oldSize = size;
        for (int i=size-1; i > 0; i--) {
            T next = pop();
            array[i] = next;
        }
        size = oldSize;
        List<T> list = new ArrayList<T>(size);
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
    public String toString() {
        Comparable<T>[] output = new Comparable[size];
        System.arraycopy(array, 0, output, 0, size);
        return Arrays.toString(output);
    }

    public static void main(String[] args) {
        FixedHeap<Integer> heap = new FixedHeap<Integer>(1, 2, 3, 4, 5, 6, 7);
        System.out.println(heap);
        for (int i=8; i<33; i++) {
            heap.push(i);
            System.out.println(heap);
        }
        System.out.println(heap.sort());
        System.out.println(heap);

        FixedHeap<String> stringHeap = new FixedHeap<String>("a", "b", "c", "d");
        System.out.println(stringHeap);
        stringHeap.push("x");
        System.out.println(stringHeap);
        System.out.println(stringHeap.sort());
        System.out.println(stringHeap);
    }
}
