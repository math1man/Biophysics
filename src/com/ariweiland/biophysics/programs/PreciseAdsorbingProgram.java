package com.ariweiland.biophysics.programs;

import acm.program.ConsoleProgram;
import com.ariweiland.biophysics.peptide.Polypeptide;
import com.ariweiland.biophysics.peptide.Residue;
import com.ariweiland.biophysics.sampler.BruteForcePreciseSurfaceSampler;
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
public class PreciseAdsorbingProgram extends ConsoleProgram {

    private final JButton fibonacci = new JButton("Fibonacci");
    private final JTextField fibN = new JTextField("7");

    private final JButton random = new JButton("Random");
    private final JTextField length = new JTextField("20");
    private final JTextField ratio = new JTextField("0.4");

    private final ButtonGroup group = new ButtonGroup();
    private final JRadioButton dim2 = new JRadioButton("2D");
    private final JRadioButton dim3 = new JRadioButton("3D");
    private int dimension = 2;

    private final JTextField surfaceEnergy = new JTextField("0.08");

    private final JButton start = new JButton("Start");
    private final JButton stop = new JButton("Stop");
    private final JButton clear = new JButton("Clear");
    private final JButton time = new JButton("Current Runtime");

    private final JTextField sequence = new JTextField("Sequence", 60);

    private long startTime = System.currentTimeMillis();
    private MyThread thread;

    @Override
    public void init() {
        setSize(800, 600);

        Residue.setInteractionScheme(-1, -1.0 / 7, 0);

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

        add(surfaceEnergy, WEST);
        surfaceEnergy.addActionListener(this);

        add(start, WEST);
        add(stop, WEST);
        add(clear, WEST);

        add(time, WEST);

        add(sequence, SOUTH);
        sequence.addActionListener(this);

        dim2.doClick();
        sequence.grabFocus();

        addActionListeners();
    }

    private double getSurfaceEnergy() {
        try {
            double value = Double.valueOf(surfaceEnergy.getText());
            if (value < 0) {
                JOptionPane.showMessageDialog(null, "The energy must be greater than or equal to zero.",
                        "Invalid Energy", JOptionPane.ERROR_MESSAGE);
            } else {
                return value;
            }
        } catch (NumberFormatException e1) {
            JOptionPane.showMessageDialog(null, "The surface energy is not a number.", "Invalid Energy", JOptionPane.ERROR_MESSAGE);
        }
        return 0;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        Object source = e.getSource();
        if (source.equals(fibonacci) || source.equals(fibN)) {
            try {
                int i = Integer.valueOf(fibN.getText());
                if (i < 0) {
                    JOptionPane.showMessageDialog(null, "The Fibonacci index must be >= 0.",
                            "Invalid Index", JOptionPane.ERROR_MESSAGE);
                } else {
                    sequence.setText(Polypeptide.fibonacci(i).toString());
                }
            } catch (NumberFormatException e1) {
                JOptionPane.showMessageDialog(null, "The Fibonacci index must be an integer value.",
                        "Invalid Index", JOptionPane.ERROR_MESSAGE);
            }
        } else if (source.equals(random) || source.equals(length) || source.equals(ratio)) {
            try {
                int l = Integer.valueOf(length.getText());
                double r = Double.valueOf(ratio.getText());
                if (l < 3) {
                    JOptionPane.showMessageDialog(null, "The minimum sequence length must be > 2.",
                            "Invalid Length", JOptionPane.ERROR_MESSAGE);
                } else if (r > 1 || r < 0) {
                    JOptionPane.showMessageDialog(null, "The H-P ratio must be between 0 and 1.",
                            "Invalid Ratio", JOptionPane.ERROR_MESSAGE);
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
        } else if (source.equals(start) || source.equals(sequence)) {
            try {
                Polypeptide polypeptide = new Polypeptide(sequence.getText());
                String[] split = surfaceEnergy.getText().split(",");
                List<Double> values = new ArrayList<>();
                for (String s : split) {
                    values.add(Double.valueOf(s));
                }
                if (polypeptide.size() < 3) {
                    JOptionPane.showMessageDialog(null, "Polypeptide '" + polypeptide + "' is invalid. The polypeptide must be at least 3 peptides long.",
                            "Invalid Polypeptide", JOptionPane.ERROR_MESSAGE);
                } else {
                    thread = new MyThread(dimension, polypeptide, values);
                    startTime = System.currentTimeMillis();
                    thread.start();
                }
            } catch (NumberFormatException e1) {
                JOptionPane.showMessageDialog(null, "The surface energy list is not all numbers.", "Invalid Energy", JOptionPane.ERROR_MESSAGE);
            }
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

    private void setFoldingEnabled(boolean isEnabled) {
        fibonacci.setEnabled(isEnabled);
        fibN.setEnabled(isEnabled);
        random.setEnabled(isEnabled);
        length.setEnabled(isEnabled);
        dim2.setEnabled(isEnabled);
        dim3.setEnabled(isEnabled);
        ratio.setEnabled(isEnabled);
        start.setEnabled(isEnabled);
        sequence.setEnabled(isEnabled);
    }

    private class MyThread extends Thread {

        private final Sampler sampler = new BruteForcePreciseSurfaceSampler(Residue.S);

        private final int dimension;
        private final Polypeptide polypeptide;
        private final List<Double> surfaceEnergies;

        private boolean running = true;

        public MyThread(int dimension, Polypeptide polypeptide, List<Double> surfaceEnergies) {
            this.dimension = dimension;
            this.polypeptide = polypeptide;
            this.surfaceEnergies = surfaceEnergies;
        }

        public void terminate() {
            running = false;
            sampler.terminate();
            setFoldingEnabled(true);
        }

        @Override
        public void run() {
            setFoldingEnabled(false);
            println(polypeptide);
            for (double se : surfaceEnergies) {
                Residue.setSurfaceInteractions(-se, -se); // attractive interactions
                println();
                println(String.format("Attractive Surface Interaction: %.3f", se));
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
                if (!running) break;
            }
            terminate();
        }
    }

    public static void main(String[] args) {
        new PreciseAdsorbingProgram().start();
    }

}
