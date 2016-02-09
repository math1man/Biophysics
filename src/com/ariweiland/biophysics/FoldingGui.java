package com.ariweiland.biophysics;

import acm.program.ConsoleProgram;
import com.ariweiland.biophysics.lattice.CheckedLattice;
import com.ariweiland.biophysics.modeler.CurrentParallelModeler;
import com.ariweiland.biophysics.modeler.CurrentSurfaceModeler;
import com.ariweiland.biophysics.modeler.Modeler;
import com.ariweiland.biophysics.peptide.Polypeptide;
import com.ariweiland.biophysics.peptide.Residue;

import javax.swing.*;
import java.awt.event.ActionEvent;

/**
 * TODO: add dimension
 *
 * @author Ari Weiland
 */
public class FoldingGui extends ConsoleProgram {

    private final JTextField sequence = new JTextField("HP_+m", 40);

    private final JCheckBox surface = new JCheckBox("Surface Fold?");
    private final ButtonGroup group = new ButtonGroup();
    private final JRadioButton hydrophobic = new JRadioButton("Hydrophobic (H)");
    private final JRadioButton hydrophylic = new JRadioButton("Hydrophilic (P)");
    private final JRadioButton positive = new JRadioButton("Positive (+)");
    private final JRadioButton negative = new JRadioButton("Negative (-)");
    private Residue surfaceType;

    private final JButton fold = new JButton("Fold");
    private final JButton stop = new JButton("Stop");
    private final JButton clear = new JButton("Clear");

    private GuiThread thread;

    @Override
    public void init() {
        setSize(800, 600);

        add(sequence, SOUTH);
        sequence.addActionListener(this);

        add(fold, SOUTH);
        add(stop, SOUTH);
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
            Modeler modeler;
            if (surface.isSelected()) {
                modeler = new CurrentSurfaceModeler(2, surfaceType);
            } else {
                modeler = new CurrentParallelModeler(2);
            }
            thread = new GuiThread(new Polypeptide(sequence.getText()), modeler);
            thread.start();
        } else if (source.equals(stop)) {
            if (thread != null) {
                thread.terminate();
            }
        } else if (source.equals(clear)) {
            getConsole().clear();
        }
    }

    private class GuiThread extends Thread {

        private final Polypeptide polypeptide;
        private final Modeler modeler;

        private GuiThread(Polypeptide polypeptide, Modeler modeler) {
            this.polypeptide = polypeptide;
            this.modeler = modeler;
        }

        public void terminate() {
            modeler.terminate();
            fold.setEnabled(true);
            sequence.setEnabled(true);
        }

        @Override
        public void run() {
            fold.setEnabled(false);
            sequence.setEnabled(false);

            println("Folding " + polypeptide + "...");
            println("Node count: " + polypeptide.size());
            println();

            long start = System.currentTimeMillis();
            CheckedLattice lattice = modeler.fold(polypeptide);
            long elapsed = System.currentTimeMillis() - start;

            for (String line : lattice.visualize()) {
                println(line);
            }

            println("Elapsed time: " + (elapsed / 1000.0) + " s");
            println("Lattice energy: " + lattice.getEnergy());
            println("Perimeter: " + lattice.getSurfaceSize() + "/" + lattice.boundingPerimeter());
            println();
            terminate();
        }
    }

    public static void main(String[] args) {
        new FoldingGui().start();
    }
}
