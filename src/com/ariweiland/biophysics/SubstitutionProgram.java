package com.ariweiland.biophysics;

import acm.program.ConsoleProgram;
import com.ariweiland.biophysics.peptide.Polypeptide;
import com.ariweiland.biophysics.peptide.Residue;
import com.ariweiland.biophysics.sampler.BruteForceSampler;
import com.ariweiland.biophysics.sampler.Sampler;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * @author Ari Weiland
 */
public class SubstitutionProgram extends ConsoleProgram {

    private final JButton fibonacci = new JButton("Fibonacci Sequence");
    private final JTextField fibN = new JTextField("7");

    private final JButton random = new JButton("Random Sequence");
    private final JTextField length = new JTextField("20");
    private final JTextField ratio = new JTextField("0.4");

    private final ButtonGroup group = new ButtonGroup();
    private final JRadioButton dim2 = new JRadioButton("2D");
    private final JRadioButton dim3 = new JRadioButton("3D");
    private int dimension = 2;

    private final JButton fold = new JButton("Start");
    private final JButton stop = new JButton("Stop");
    private final JButton clear = new JButton("Clear");
    private final JButton time = new JButton("Current Runtime");

    private final JTextField sequence = new JTextField("Sequence", 60);

    private long startTime = System.currentTimeMillis();
    private MyThread thread;

    @Override
    public void init() {
        setSize(800, 600);

        Residue.setInteractionScheme(-1, 0, 0);

        add(fibonacci, WEST);
        add(fibN, WEST);
        fibN.addActionListener(this);

        add(random, WEST);
        add(length, WEST);
        length.addActionListener(this);
        add(ratio, WEST);
        ratio.addActionListener(this);

        group.add(dim2);
        add(dim2, WEST);
        dim2.addActionListener(this);
        group.add(dim3);
        add(dim3, WEST);
        dim3.addActionListener(this);

        add(fold, WEST);
        add(stop, WEST);
        add(clear, WEST);

        add(time, WEST);

        add(sequence, SOUTH);
        sequence.addActionListener(this);

        dim2.doClick();
        sequence.grabFocus();

        addActionListeners();
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        Object source = e.getSource();
        if (source.equals(fibonacci) || source.equals(fibN)) {
            try {
                int i = Integer.valueOf(fibN.getText());
                if (i < 0) {
                    JOptionPane.showMessageDialog(null, "The Fibonacci index must be >= 0.", "Invalid Index", JOptionPane.ERROR_MESSAGE);
                } else {
                    sequence.setText(Polypeptide.fibonacci(i).toString());
                }
            } catch (NumberFormatException e1) {
                JOptionPane.showMessageDialog(null, "The Fibonacci index must be an integer value.", "Invalid Index", JOptionPane.ERROR_MESSAGE);
            }
        } else if (source.equals(random) || source.equals(length) || source.equals(ratio)) {
            try {
                int l = Integer.valueOf(length.getText());
                double r = Double.valueOf(ratio.getText());
                if (l < 3) {
                    JOptionPane.showMessageDialog(null, "The minimum sequence length must be > 2.", "Invalid Length", JOptionPane.ERROR_MESSAGE);
                } else if (r > 1 || r < 0) {
                    JOptionPane.showMessageDialog(null, "The H-P ratio must be between 0 and 1.", "Invalid Ratio", JOptionPane.ERROR_MESSAGE);
                } else {
                    sequence.setText(Polypeptide.random(l, r).toString());
                }
            } catch (NumberFormatException e1) {
                JOptionPane.showMessageDialog(null, "One of the random sequence parameters is not a number.", "Invalid Parameters", JOptionPane.ERROR_MESSAGE);
            }
        } else if (source.equals(dim2)) {
            dimension = 2;
        } else if (source.equals(dim3)) {
            dimension = 3;
        } else if (source.equals(fold)) {
            thread = new MyThread(dimension, new Polypeptide(sequence.getText()));
            startTime = System.currentTimeMillis();
            thread.start();
        } else if (source.equals(stop)) {
            if (thread != null) {
                thread.terminate();
            }
        } else if (source.equals(clear)) {
            getConsole().clear();
        } else if (source.equals(time)) {
            long elapsed = System.currentTimeMillis() - startTime;
            println(String.format("\n>>> Current Running Time: %.2f minutes\n", elapsed / 60000.0));
        }
    }

    private class MyThread extends Thread {

        private final int dimension;
        private final Polypeptide original;
        private final BruteForceSampler sampler = new BruteForceSampler();

        private boolean running = true;

        public MyThread(int dimension, Polypeptide original) {
            this.dimension = dimension;
            this.original = original;
        }

        public void terminate() {
            running = false;
            sampler.terminate();
            fold.setEnabled(true);
        }

        private Polypeptide changePeptide(Polypeptide p, int i) {
            char[] sequence = p.toString().replaceAll("[()-]", "").toCharArray();
            if (p.get(i).residue == Residue.P) {
                sequence[i] = 'H';
            } else {
                sequence[i] = 'P';
            }
            return new Polypeptide(String.valueOf(sequence));
        }

        @Override
        public void run() {
            fold.setEnabled(false);
            for (int i=-1; i<original.size() && running; i++) {
                Polypeptide polypeptide;
                if (i < 0) {
                    polypeptide = original;
                    println("Unsubstituted:");
                } else {
                    polypeptide = changePeptide(original, i);
                    println("Substituting position " + i + ":");
                }
                println(polypeptide);
                long start = System.currentTimeMillis();
                Map<Double, Double> density = sampler.normalize(sampler.getDensity(dimension, polypeptide));
                long elapsed = System.currentTimeMillis() - start;
                println("Elapsed time: " + (elapsed / 1000.0) + " s");
                println("Bins\tCounts");
                List<Double> keys = new ArrayList<>(density.keySet());
                Collections.sort(keys);
                for (double d : keys) {
                    println(d + " \t" + density.get(d));
                }
                println();
                println(Sampler.asMathematicaCode(density));
                println();
            }
            terminate();
        }
    }

    public static void main(String[] args) {
        new SubstitutionProgram().start();
    }
}
