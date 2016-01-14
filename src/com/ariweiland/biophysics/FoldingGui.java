package com.ariweiland.biophysics;

import acm.program.ConsoleProgram;
import com.ariweiland.biophysics.lattice.Lattice;
import com.ariweiland.biophysics.modeler.CurrentParallelModeler;
import com.ariweiland.biophysics.modeler.CurrentSurfaceModeler;
import com.ariweiland.biophysics.modeler.Modeler;
import com.ariweiland.biophysics.peptide.Polypeptide;
import com.ariweiland.biophysics.peptide.Residue;

import javax.swing.*;
import java.awt.event.ActionEvent;

/**
 * @author Ari Weiland
 */
public class FoldingGui extends ConsoleProgram {

    private final JTextField sequence = new JTextField("HP +m", 30);

    private final JCheckBox surface = new JCheckBox("Surface Fold?");
    private final ButtonGroup group = new ButtonGroup();
    private final JRadioButton hydrophobic = new JRadioButton("Hydrophobic (H)");
    private final JRadioButton hydrophylic = new JRadioButton("Hydrophilic (P)");
    private final JRadioButton positive = new JRadioButton("Positive (+)");
    private final JRadioButton negative = new JRadioButton("Negative (-)");
    private Residue surfaceType;

    private final JButton fold = new JButton("Fold");
    private final JButton clear = new JButton("Clear");

    @Override
    public void init() {
        setSize(800, 600);

        add(sequence, SOUTH);
        sequence.addActionListener(this);

        add(fold, SOUTH);
        add(clear, SOUTH);

        add(surface, EAST);
        surface.addActionListener(this);

        group.add(hydrophobic);
        add(hydrophobic, EAST);
        hydrophobic.addActionListener(this);
        group.add(hydrophylic);
        add(hydrophylic, EAST);
        hydrophylic.addActionListener(this);
        group.add(positive);
        add(positive, EAST);
        positive.addActionListener(this);
        group.add(negative);
        add(negative, EAST);
        negative.addActionListener(this);

        hydrophobic.doClick();
        sequence.grabFocus();

        addActionListeners();
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        Object source = e.getSource();
        if (source.equals(hydrophylic)) {
            surfaceType = Residue.P;
        } else if (source.equals(hydrophobic)) {
            surfaceType = Residue.H;
        } else if (source.equals(positive)) {
            surfaceType = Residue.POS;
        } else if (source.equals(negative)) {
            surfaceType = Residue.NEG;
        } else if (source.equals(sequence) || source.equals(fold)) {
            String temp = sequence.getText();
            clear();
            foldPolypeptide(temp);
        } else if (source.equals(clear)) {
            clear();
        }
    }

    public void foldPolypeptide(String sequence) {
        Polypeptide polypeptide = new Polypeptide(sequence);
        Modeler modeler;
        if (surface.isSelected()) {
            modeler = new CurrentSurfaceModeler(surfaceType);
        } else {
            modeler = new CurrentParallelModeler();
        }

        println(polypeptide);
        println("Node count: " + polypeptide.size());
        println();

        long start = System.currentTimeMillis();
        Lattice lattice = modeler.fold(polypeptide);
        long elapsed = System.currentTimeMillis() - start;

        for (String line : lattice.visualize()) {
            println(line);
        }

        println("Elapsed time: " + (elapsed / 1000.0) + " s");
        println("Lattice energy: " + lattice.getEnergy());
        println("Perimeter: " + lattice.getPerimeter() + "/" + lattice.boundingPerimeter());
    }

    public void clear() {
        getConsole().clear();
    }

    public static void main(String[] args) {
        new FoldingGui().start();
    }
}
