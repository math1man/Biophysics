package com.ariweiland.biophysics;

import acm.program.ConsoleProgram;
import com.ariweiland.biophysics.peptide.Polypeptide;
import com.ariweiland.biophysics.peptide.Residue;
import com.ariweiland.biophysics.sampler.BruteForceSampler;
import com.ariweiland.biophysics.sampler.DefaultWangLandauSampler;
import com.ariweiland.biophysics.sampler.NaiveSampler;
import com.ariweiland.biophysics.sampler.WangLandauSampler;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * @author Ari Weiland
 */
public class SamplerProgram extends ConsoleProgram {

    private final JTextField ratio = new JTextField("H-P Ratio");
    private final JTextField minLength = new JTextField("Min Sequence Length");
    private final JTextField maxLength = new JTextField("Max Sequence Length");

    private final ButtonGroup group = new ButtonGroup();
    private final JRadioButton dim2 = new JRadioButton("2D");
    private final JRadioButton dim3 = new JRadioButton("3D");
    private int dimension = 2;

    private final JButton fold = new JButton("Start");
    private final JButton stop = new JButton("Stop");
    private final JButton clear = new JButton("Clear");

    private final JCheckBox repeat = new JCheckBox("Repeat?");

    private final JButton time = new JButton("Current Runtime");

    private long startTime = System.currentTimeMillis();
    private MyThread thread;

    @Override
    public void init() {
        setSize(800, 600);

        add(ratio, WEST);
        ratio.addActionListener(this);
        add(minLength, WEST);
        minLength.addActionListener(this);
        add(maxLength, WEST);
        maxLength.addActionListener(this);

        group.add(dim2);
        add(dim2, WEST);
        dim2.addActionListener(this);
        group.add(dim3);
        add(dim3, WEST);
        dim3.addActionListener(this);

        add(fold, WEST);
        add(stop, WEST);
        add(clear, WEST);

        add(repeat, WEST);
        repeat.addActionListener(this);

        add(time, WEST);

        dim2.doClick();
        ratio.grabFocus();
        
        addActionListeners();
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        Object source = e.getSource();
        if (source.equals(dim2)) {
            dimension = 2;
        } else if (source.equals(dim3)) {
            dimension = 3;
        } else if (source.equals(fold)) {
            thread = new MyThread(dimension, repeat.isSelected(), ratio.getText(), minLength.getText(), maxLength.getText());
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
        private boolean repeat;
        private double ratio = Double.NaN;
        private int minLength = 0;
        private int maxLength = 0;
        private final BruteForceSampler bfs = new BruteForceSampler();
        private final NaiveSampler ns = new NaiveSampler();
        private final WangLandauSampler wls = new DefaultWangLandauSampler();

        private boolean running = true;

        public MyThread(int dimension, boolean repeat, String ratioString, String minString, String maxString) {
            this.dimension = dimension;
            this.repeat = repeat;
            try {
                ratio = Double.valueOf(ratioString);
                if (ratio > 1 || ratio < 0) {
                    JOptionPane.showMessageDialog(null, "The H-P ratio must be between 0 and 1.", "Invalid Ratio", JOptionPane.ERROR_MESSAGE);
                    running = false;
                    return;
                }
            } catch (NumberFormatException e1) {
                JOptionPane.showMessageDialog(null, "The H-P ratio must be a decimal value.", "Invalid Ratio", JOptionPane.ERROR_MESSAGE);
                running = false;
                return;
            }
            try {
                minLength = Integer.valueOf(minString);
                if (minLength < 3) {
                    JOptionPane.showMessageDialog(null, "The minimum sequence length must be > 2.", "Invalid Length", JOptionPane.ERROR_MESSAGE);
                    running = false;
                    return;
                }
            } catch (NumberFormatException e1) {
                JOptionPane.showMessageDialog(null, "The minimum sequence length must be an integer value.", "Invalid Length", JOptionPane.ERROR_MESSAGE);
                running = false;
                return;
            }
            try {
                maxLength = Integer.valueOf(maxString);
                if (maxLength < 3) {
                    JOptionPane.showMessageDialog(null, "The maximum sequence length must be > 2.", "Invalid Length", JOptionPane.ERROR_MESSAGE);
                    running = false;
                    return;
                }
            } catch (NumberFormatException e1) {
                JOptionPane.showMessageDialog(null, "The maximum sequence length must be an integer value.", "Invalid Length", JOptionPane.ERROR_MESSAGE);
                running = false;
                return;
            }
            if (minLength > maxLength) {
                JOptionPane.showMessageDialog(null, "The maximum sequence length must be not be less than the minimum sequence",
                        "Invalid Length", JOptionPane.ERROR_MESSAGE);
                running = false;
            }
        }

        public void terminate() {
            repeat = false;
            running = false;
            bfs.terminate();
            ns.terminate();
            wls.terminate();
            fold.setEnabled(true);
        }

        @Override
        public void run() {
            while (running) {
                fold.setEnabled(false);
                println("H-P Ratio:   " + ratio);
                println("Min Seq Len: " + minLength);
                println("Max Seq Len: " + maxLength);
                println();
                for (int i=minLength; i<=maxLength && running; i++) {
                    List<Residue> residues = new ArrayList<>();
                    for (int j=0; j<i; j++) {
                        if (j < (i * ratio)) {
                            residues.add(Residue.P);
                        } else {
                            residues.add(Residue.H);
                        }
                    }
                    Collections.shuffle(residues);
                    Polypeptide polypeptide = new Polypeptide(residues);
                    println("Folding " + polypeptide + "...");
                    println("Node count: " + polypeptide.size());
                    println();

                    println("Brute Force:");
                    long start = System.currentTimeMillis();
                    Map<Double, Double> density = bfs.normalize(bfs.getDensity(dimension, polypeptide));
                    long elapsed = System.currentTimeMillis() - start;
                    println("Elapsed time: " + (elapsed / 1000.0) + " s");
                    println("Energy\tDensities");
                    List<Double> keys = new ArrayList<>(density.keySet());
                    Collections.sort(keys);
                    for (double d : keys) {
                        println(d + " \t" + density.get(d));
                    }
                    println();

                    if (running) {
                        println("Naive:");
                        start = System.currentTimeMillis();
                        density = ns.normalize(ns.getDensity(dimension, polypeptide));
                        elapsed = System.currentTimeMillis() - start;
                        println("Elapsed time: " + (elapsed / 1000.0) + " s");
                        println("Energy\tDensities");
                        keys = new ArrayList<>(density.keySet());
                        Collections.sort(keys);
                        for (double d : keys) {
                            println(d + " \t" + density.get(d));
                        }
                        println();
                    }

                    if (running) {
                        println("Wang-Landau:");
                        start = System.currentTimeMillis();
                        density = wls.normalize(wls.getDensity(dimension, polypeptide));
                        elapsed = System.currentTimeMillis() - start;
                        println("Elapsed time: " + (elapsed / 1000.0) + " s");
                        println("Energy\tDensities");
                        keys = new ArrayList<>(density.keySet());
                        Collections.sort(keys);
                        for (double d : keys) {
                            println(d + " \t" + density.get(d));
                        }
                        println();
                    }
                    println("=============================");
                    println();
                }
                running = repeat;
            }
            terminate();
        }
    }

    public static void main(String[] args) {
        new SamplerProgram().start();
    }
}
