package net.calm.embryogen.simulate;

import ij.IJ;
import ij.ImageStack;
import net.calm.embryogen.image_simulator.Nucleus;
import net.calm.embryogen.intensity.IntensThread;
import net.calm.embryogen.params.SimParams;

import java.util.Random;

public class Simulation {
    private final SimParams params;
    private final double eps;
    private final double gam;
    private final double Lx;
    private final double Ly;
    private final double Lz;
    private final Random r;
    private final double PI = Math.PI;
    private int count_BM;

    public Simulation(SimParams params, double eps, double gam, double[] dims) {
        this.params = params;
        this.eps = eps;
        this.gam = gam;
        this.Lx = dims[0];
        this.Ly = dims[1];
        this.Lz = dims[2];
        this.r = new Random();
        this.count_BM = 0;
    }

    public void simulation(Nucleus[] a, ImageStack output) {
        if (params.isCluster()) initial_relaxation(a);
        Tanh_blob(a, output);
    }

    void initial_relaxation(Nucleus[] a) {
        int i, imax;
        double Dxc, Dyc, Dzc;

        Dxc = params.getDx();
        Dyc = params.getDy();
        Dzc = params.getDz();

        imax = (int) Math.round(500.0 / params.getDT());
        params.setDx(0.0);
        params.setDy(0.0);
        params.setDz(0.0);

        for (i = 0; i < imax; i++) {
            Euler_Method(a);
        }

        params.setDx(Dxc);
        params.setDy(Dyc);
        params.setDz(Dzc);

    }

    void Tanh_blob(Nucleus[] a, ImageStack output) {
        int nbCPUs = Runtime.getRuntime().availableProcessors();
        IntensThread[] intensThreads = new IntensThread[nbCPUs];
        for (int thread = 0; thread < nbCPUs; thread++) {
            intensThreads[thread] = new IntensThread(thread, nbCPUs, a, output, params);
            intensThreads[thread].start();
        }
        try {
            for (int thread = 0; thread < nbCPUs; thread++) {
                intensThreads[thread].join();
            }
        } catch (InterruptedException ie) {
            IJ.error("A thread was interrupted during output generation.");
        }

    }

    void Euler_Method(Nucleus[] a) {
        Nucleus[] temp = new Nucleus[a.length];

        double xi;
        double xij, yij, zij, rij, Fijx, Fijy, Fijz;
        int i, j;

        for (i = 0; i < a.length; i++) {
            temp[i] = a[i];
        }

        for (i = 0; i < a.length; i++) {

            Fijx = 0.0;
            Fijy = 0.0;
            Fijz = 0.0;

            for (j = 0; j < a.length; j++) {
                if (j != i) {
                    xij = a[j].getX() - a[i].getX();
                    yij = a[j].getY() - a[i].getY();
                    zij = a[j].getZ() - a[i].getZ();

                    rij = Math.sqrt(xij * xij + yij * yij + zij * zij);
                    Fijx += eps * (-3.0 * Math.pow(gam, 3.0) / Math.pow(rij, 5.0) + gam / Math.pow(rij, 3.0)) * (-xij);
                    Fijy += eps * (-3.0 * Math.pow(gam, 3.0) / Math.pow(rij, 5.0) + gam / Math.pow(rij, 3.0)) * (-yij);
                    Fijz += eps * (-3.0 * Math.pow(gam, 3.0) / Math.pow(rij, 5.0) + gam / Math.pow(rij, 3.0)) * (-zij);
                }

            }

            xi = Box_Muller_Method(0.0, 1.0);
            temp[i].setX(a[i].getX() - Fijx * params.getDT() + Math.sqrt(2.0 * params.getDx()) * Math.sqrt(params.getDT()) * xi);

            if (temp[i].getX() < 0) {
                temp[i].setX(0.0);
            }
            if (temp[i].getX() > Lx) {
                temp[i].setX(Lx);
            }

            xi = Box_Muller_Method(0.0, 1.0);
            temp[i].setY(a[i].getY() - Fijy * params.getDT() + Math.sqrt(2.0 * params.getDy()) * Math.sqrt(params.getDT()) * xi);

            if (temp[i].getY() < 0) {
                temp[i].setY(0.0);
            }
            if (temp[i].getY() > Ly) {
                temp[i].setY(Ly);
            }

            xi = Box_Muller_Method(0.0, 1.0);
            temp[i].setZ(a[i].getZ() - Fijz * params.getDT() + Math.sqrt(2.0 * params.getDz()) * Math.sqrt(params.getDT()) * xi);

            if (temp[i].getZ() < 0) {
                temp[i].setZ(0.0);
            }
            if (temp[i].getZ() > Lz) {
                temp[i].setZ(Lz);
            }

        }

        for (i = 0; i < a.length; i++) {
            a[i] = temp[i];
        }
    }

    double Box_Muller_Method(double m, double sd) {

        double r1, r2;
        double nrand;

        r1 = r.nextDouble();
        r2 = r.nextDouble();

        do {

            if (count_BM == 0) {
                count_BM = 1;
                nrand = sd * Math.sqrt(-2 * Math.log(r1)) * Math.cos(2 * PI * r2) + m;
            } else {
                count_BM = 0;
                nrand = sd * Math.sqrt(-2 * Math.log(r1)) * Math.sin(2 * PI * r2) + m;
            }

        } while (Double.isInfinite(Math.sqrt(-2 * Math.log(r1))));

        return nrand;

    }

}
