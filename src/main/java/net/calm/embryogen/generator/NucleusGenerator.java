package net.calm.embryogen.generator;

import net.calm.embryogen.image_simulator.Nucleus;
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

    public void initialiseNuclei(Nucleus[] a) {
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Initialising nuclei..."));
        //setting the initial positions of particle
        for (int i = 0; i < nCells; i++) {
            a[i] = new Nucleus();
            a[i].setTheta_x(0.25 * PI * r.nextFloat());
            a[i].setTheta_y(0.25 * PI * r.nextFloat());
            a[i].setTheta_z(0.5 * PI * r.nextFloat());
            if (params.isCluster()) {
                a[i].setX(0.4 * Lx + 0.2 * Lx * r.nextDouble());
                a[i].setY(0.4 * Ly + 0.2 * Ly * r.nextDouble());
                a[i].setZ(0.4 * Lz + 0.2 * Lz * r.nextDouble());
            } else {
                a[i].setX(Lx * r.nextDouble());
                a[i].setY(Ly * r.nextDouble());
                a[i].setZ(Lz * r.nextDouble());
            }
            a[i].setxLength(params.getA() + 0.2 * params.getA() * r.nextGaussian());
            a[i].setyLength(params.getB() + 0.2 * params.getB() * r.nextGaussian());
            a[i].setzLength(params.getC() + 0.2 * params.getC() * r.nextGaussian());
        }
    }
}
