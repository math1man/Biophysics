package com.ariweiland.biophysics.programs;

import acm.program.ConsoleProgram;
import com.ariweiland.biophysics.peptide.Polypeptide;
import com.ariweiland.biophysics.peptide.Residue;
import com.ariweiland.biophysics.sampler.BruteForceSurfaceSampler;
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
public class AdsorbingProgram extends ConsoleProgram {


    private final JButton fibonacci = new JButton("Fibonacci");
    private final JTextField fibN = new JTextField("7");

    private final JButton random = new JButton("Random");
    private final JTextField length = new JTextField("20");
    private final JTextField ratio = new JTextField("0.4");

    private final ButtonGroup group = new ButtonGroup();
    private final JRadioButton dim2 = new JRadioButton("2D");
    private final JRadioButton dim3 = new JRadioButton("3D");
    private int dimension = 2;

    private final JTextField surfaceEnergyRange = new JTextField("0:0.05:1");

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

        Residue.setInteractionScheme(-1, -1.0/7, 0);

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

        add(surfaceEnergyRange, WEST);
        surfaceEnergyRange.addActionListener(this);

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

    private double[] getSurfaceEnergyRange() {
        String[] split = surfaceEnergyRange.getText().split(":");
        if (split.length == 3) {
            try {
                double min = Integer.valueOf(split[0]);
                double increment = Double.valueOf(split[1]);
                double max = Double.valueOf(split[2]);
                if (min > max) {
                    JOptionPane.showMessageDialog(null, "The minimum energy must be less than or equal to the max energy.",
                            "Invalid Energies", JOptionPane.ERROR_MESSAGE);
                } else {
                    return new double[]{min, increment, max};
                }
            } catch (NumberFormatException e1) {
                JOptionPane.showMessageDialog(null, "One of the random sequence parameters is not a number.", "Invalid Parameters", JOptionPane.ERROR_MESSAGE);
            }
        } else {
            JOptionPane.showMessageDialog(null, "The surface energy range must be of the form '#1:#2:#3', " +
                    "where #1 is the minimum value, #2 is the increment, and #3 is the maximum value.",
                    "Invalid Energy Range", JOptionPane.ERROR_MESSAGE);
        }
        return new double[]{0,0,0};
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
            String[] split = surfaceEnergyRange.getText().split(":");
            if (split.length == 3) {
                try {
                    Polypeptide polypeptide = new Polypeptide(sequence.getText());
                    double min = Double.valueOf(split[0]);
                    double increment = Double.valueOf(split[1]);
                    double max = Double.valueOf(split[2]);
                    if (polypeptide.size() < 3) {
                        JOptionPane.showMessageDialog(null, "Polypeptide '" + polypeptide + "' is invalid. The polypeptide must be at least 3 peptides long.",
                                "Invalid Polypeptide", JOptionPane.ERROR_MESSAGE);
                    } else if (min > max) {
                        JOptionPane.showMessageDialog(null, "The minimum energy must be less than or equal to the max energy.",
                                "Invalid Energies", JOptionPane.ERROR_MESSAGE);
                    } else if (increment <= 0) {
                        JOptionPane.showMessageDialog(null, "The energy increment must be greater than 0.",
                                "Invalid Increment", JOptionPane.ERROR_MESSAGE);
                    } else {
                        thread = new MyThread(dimension, polypeptide, min, increment, max);
                        startTime = System.currentTimeMillis();
                        thread.start();
                    }
                } catch (NumberFormatException e1) {
                    JOptionPane.showMessageDialog(null, "One of the random sequence parameters is not a number.", "Invalid Parameters", JOptionPane.ERROR_MESSAGE);
                }
            } else {
                JOptionPane.showMessageDialog(null, "The surface energy range must be of the form '#1:#2:#3', " +
                                "where #1 is the minimum value, #2 is the increment, and #3 is the maximum value.",
                        "Invalid Energy Range", JOptionPane.ERROR_MESSAGE);
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

        private final Sampler sampler = new BruteForceSurfaceSampler(Residue.S);

        private final int dimension;
        private final Polypeptide polypeptide;
        private final double minSurfaceEnergy;
        private final double surfaceEnergyIncrement;
        private final double maxSurfaceEnergy;

        private boolean running = true;

        public MyThread(int dimension, Polypeptide polypeptide, double minSurfaceEnergy, double surfaceEnergyIncrement, double maxSurfaceEnergy) {
            this.dimension = dimension;
            this.polypeptide = polypeptide;
            this.minSurfaceEnergy = minSurfaceEnergy;
            this.surfaceEnergyIncrement = surfaceEnergyIncrement;
            this.maxSurfaceEnergy = maxSurfaceEnergy;
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
            for (double i = minSurfaceEnergy; i < maxSurfaceEnergy + surfaceEnergyIncrement && running; i += surfaceEnergyIncrement) {
                Residue.setSurfaceInteractions(-i, -i); // attractive interactions
                println();
                println(String.format("Attractive Surface Interaction: %.3f", i));
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
            }
            terminate();
        }
    }

    public static void main(String[] args) {
        new AdsorbingProgram().start();
    }

}
