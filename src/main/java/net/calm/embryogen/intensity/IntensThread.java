package net.calm.embryogen.intensity;

import ij.ImageStack;
import net.calm.embryogen.image_simulator.Nucleus;
import net.calm.embryogen.params.SimParams;

public class IntensThread extends Thread {

    private final Nucleus[] a;
    private final ImageStack output;
    private final int thread;
    private final int nThreads;
    private final SimParams params;

    public IntensThread(int thread, int nThreads, Nucleus[] a, ImageStack output, SimParams params) {
        this.thread = thread;
        this.nThreads = nThreads;
        this.a = a;
        this.output = output;
        this.params = params;
    }

    public void run() {
        for (int k = thread; k < output.getSize(); k += nThreads) {
            double varIntens = params.getNucMaxIntens() - (params.getNucMaxIntens() - params.getNucMinIntens()) * k / output.size();
            for (int j = 0; j < output.getHeight(); j++) {
                for (int i = 0; i < output.getWidth(); i++) {
                    double xi = params.getSimSizeX() * (2.0 * i + 1) / 2.0;
                    double yi = params.getSimSizeY() * (2.0 * j + 1) / 2.0;
                    double zi = params.getSimSizeZ() * (2.0 * k + 1) / 2.0;
                    double intensity = params.getIback();
                    double intmax = 0.0;
                    for (int l = 0; l < a.length; l++) {
                        double xd = (xi - a[l].getX());
                        double yd = (yi - a[l].getY());
                        double zd = (zi - a[l].getZ());

                        //z-axis rotation
                        double xr = Math.cos(a[l].getTheta_z()) * xd - Math.sin(a[l].getTheta_z()) * yd;
                        double yr = Math.sin(a[l].getTheta_z()) * xd + Math.cos(a[l].getTheta_z()) * yd;
                        double zr = zd;

                        //y-axis rotation
                        double zr2 = Math.cos(a[l].getTheta_y()) * zr - Math.sin(a[l].getTheta_y()) * xr;
                        double xr2 = Math.sin(a[l].getTheta_y()) * zr + Math.cos(a[l].getTheta_y()) * xr;
                        double yr2 = yr;

                        //x-axis rotation
                        yr = Math.cos(a[l].getTheta_x()) * yr2 - Math.sin(a[l].getTheta_x()) * zr2;
                        zr = Math.sin(a[l].getTheta_x()) * yr2 + Math.cos(a[l].getTheta_x()) * zr2;
                        xr = xr2;

                        double temp = 1 - (xr * xr / (params.getA() * params.getA()) + yr * yr / (params.getB() * params.getB()) + zr * zr / (params.getC() * params.getC()));
                        temp = (Math.tanh(params.getSigma() * temp) + 1.0) / 2.0;
                        temp = (varIntens / Math.abs((Math.tanh(params.getSigma()) + 1.0) / 2.0)) * temp;

                        if (temp > intmax) {
                            intmax = temp;
                        }
                    }
                    intensity += intmax;
                    output.setVoxel(i, j, k, intensity);
                }
            }
        }
    }
}
