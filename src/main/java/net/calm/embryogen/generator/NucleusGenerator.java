package net.calm.embryogen.generator;

import net.calm.embryogen.image_simulator.Nucleus;
import net.calm.embryogen.image_simulator.NucleusGroup;
import net.calm.embryogen.params.SimParams;
import net.calm.iaclasslibrary.TimeAndDate.TimeAndDate;

import java.util.Random;

public class NucleusGenerator {
    private final int nCells;
    private final double PI = Math.PI;
    private final Random r;
    private final double Lx;
    private final double Ly;
    private final double Lz;
    private final SimParams params;

    public NucleusGenerator(int nCells, double[] dims, SimParams params) {
        this.nCells = nCells;
        this.Lx = dims[0];
        this.Ly = dims[1];
        this.Lz = dims[2];
        this.params = params;
        r = new Random();
    }

    public void initialiseNuclei(NucleusGroup[] a) {
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Initialising nuclei..."));
        //setting the initial positions of particle
        for (int i = 0; i < nCells; i++) {
            int nChildren = r.nextDouble() < 0.01 ? (int) Math.round(Math.abs(10.0 * r.nextDouble())) : 0;
            a[i] = new NucleusGroup(
                    new Nucleus(
                            new double[]{Lx / 2.0, Ly / 2.0, Lz / 2.0},
                            new double[]{params.getA(), params.getB(), params.getC()},
                            r, params.getPositionAdjustment()),
                    nChildren);
        }
    }
}
