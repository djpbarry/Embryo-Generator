/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package image_simulator;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.GaussianBlur3D;
import ij.plugin.ImageCalculator;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.util.Random;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.image3d.regionGrowing.Watershed3D;
import org.apache.commons.math3.distribution.GammaDistribution;

/**
 *
 * @author David Barry <david.barry at crick dot ac dot uk>
 */
public class Image_Simulator {

    double Dx = 0.0, Dy = 0.0, Dz = 0.0; //Diffusion coefficient in x, y, and z directions
    double DT = 0.01; //time grid for SDE

    double px, py, pz; //pixel size in x, y, and z axis

    double framerate = 1.0; //frame rate of images

    double Lx, Ly, Lz; //domain size

    double I0, sigma = 2.0, Iback = 1000, A = 10.0, B = 5.0, C = 5.0;
    double eps = 0.1, gam = 20.0;

    int maxframe;

    float PI = 3.14159265358979323846f;

    private Random r = new Random();

    static int count_BM = 0;

    double SNR = 2.0;

    int N = 10 + r.nextInt(25);

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        new Image_Simulator();
        System.exit(0);
    }

    public Image_Simulator() {
        Nucleus[] a = new Nucleus[N];
        int tmax;
        int i;

        tmax = (int) Math.round(1.5 / DT);

        double Im = 20000; //mean intensity in nucleus

        I0 = Im - Iback;

        //voxel size
        px = 0.0962000 * 10;
        py = 0.0962000 * 10;
        pz = 0.8398000 * 10;

        //imaging domain size
        Lx = 150; // microns
        Ly = 150; //microns
        Lz = 200.0; //microns

        //number of pixels in each direction
        int nx = (int) Math.round(Lx / px);
        int ny = (int) Math.round(Ly / py);
        int nz = (int) Math.round(Lz / pz);

        ImageStack output = new ImageStack(nx, ny);
        for (int n = 1; n <= nz; n++) {
            output.addSlice(new FloatProcessor(nx, ny));
        }

        maxframe = (int) Math.round(tmax * DT * framerate);

        //setting the initial positions of particle
        for (i = 0; i < N; i++) {
            a[i] = new Nucleus();
            a[i].setTheta_x(0.25 * PI * r.nextFloat());
            a[i].setTheta_y(0.25 * PI * r.nextFloat());
            a[i].setTheta_z(0.5 * PI * r.nextFloat());

//            a[i].setX(Lx * r.nextFloat());
//            a[i].setY(Ly * r.nextFloat());
//            a[i].setZ(Lz * r.nextFloat());
            a[i].setX(0.4 * Lx + 0.2 * Lx * r.nextDouble());
            a[i].setY(0.4 * Ly + 0.2 * Ly * r.nextDouble());
            a[i].setZ(0.4 * Lz + 0.2 * Lz * r.nextDouble());

        }

        simulation(a, tmax, output);

        ImageStack pointOutput = new ImageStack(nx, ny);
        for (int n = 1; n <= nz; n++) {
            ByteProcessor slice = new ByteProcessor(nx, ny);
            slice.setValue(0.0);
            slice.fill();
            pointOutput.addSlice(slice);
        }
        for (i = 0; i < N; i++) {
            int x = (int) Math.round(a[i].getX() / px);
            int y = (int) Math.round(a[i].getY() / py);
            int z = (int) Math.round(a[i].getZ() / pz);
            pointOutput.setVoxel(x, y, z, 255);
        }
        ImageHandler intImage = ImageInt.wrap(pointOutput);
        intImage.setScale(px, pz, "microns");
        ImageHandler edt = EDT.run(intImage, 1, true, -1);
        ImageHandler edtDup = edt.duplicate();
        float edtMaxValue = (float) edtDup.getMax();
        edtDup.invert();
        edtDup.addValue(edtMaxValue + 1.0f);
        Watershed3D watershed = new Watershed3D(edtDup, intImage, 90.0, 1);
        ImageInt lineImage = watershed.getDamImage();
        ImageStack thresholdedEDT = new ImageStack(nx, ny);
        for (i = 1; i <= nz; i++) {
            ImageProcessor slice = edtDup.getImageStack().getProcessor(i);
            slice.setThreshold(0.0, 90.0, ImageProcessor.NO_LUT_UPDATE);
            slice = slice.createMask();
            assert (slice instanceof ByteProcessor);
            ((ByteProcessor) slice).outline();
            slice.invert();
            thresholdedEDT.addSlice(slice);
        }
        ImagePlus voronoiPlusMaskOutline = (new ImageCalculator()).run("Add create stack", new ImagePlus("", thresholdedEDT), lineImage.getImagePlus());
        ImageStack voronoiPlusMaskOutlineStack = new ImageStack(nx, ny);
        for (i = 1; i <= nz; i++) {
            ImageProcessor slice = voronoiPlusMaskOutline.getImageStack().getProcessor(i).convertToFloatProcessor();
            slice.multiply(1.0 / 255.0);
//            slice.add(Iback);
            voronoiPlusMaskOutlineStack.addSlice(slice);
        }
        addNoise(output);
        ImagePlus voronoiPlusEDT = (new ImageCalculator()).run("Multiply create stack", new ImagePlus("", voronoiPlusMaskOutlineStack), edt.getImagePlus());
        for (i = 1; i <= nz; i++) {
            ImageProcessor slice = voronoiPlusEDT.getImageStack().getProcessor(i);
            slice.add(1.0);
//            slice.multiply(100.0);
        }
        GaussianBlur3D.blur(voronoiPlusEDT, 2.5 / px, 2.5 / py, 2.5 / pz);
        addNoise(voronoiPlusEDT.getImageStack());
        IJ.saveAs(voronoiPlusEDT, "TIF", "D:/debugging/SimImages/voronoiPlusMaskOutline_stack.tif");
//        IJ.saveAs(new ImagePlus("", thresholdedEDT), "TIF", "D:/debugging/SimImages/thresholded_edt_stack.tif");
//        IJ.saveAs(edt.getImagePlus(), "TIF", "D:/debugging/SimImages/edt_stack.tif");
//        IJ.saveAs(lineImage.getImagePlus(), "TIF", "D:/debugging/SimImages/voronoi_lines_stack.tif");
//        IJ.saveAs(new ImagePlus("", pointOutput), "TIF", "D:/debugging/SimImages/centroid_stack.tif");
        IJ.saveAs(new ImagePlus("", output), "TIF", "D:/debugging/SimImages/sim_stack.tif");
    }

    void simulation(Nucleus[] a, int tmax, ImageStack output) {
        initial_relaxation(a);

        Tanh_blob(a, output);
    }

    void initial_relaxation(Nucleus[] a) {
        int i, imax;
        double Dxc, Dyc, Dzc;

        Dxc = Dx;
        Dyc = Dy;
        Dzc = Dz;

        imax = (int) Math.round(100.0 / DT);
        Dx = 0.0;
        Dy = 0.0;
        Dz = 0.0;

        for (i = 0; i < imax; i++) {
            Euler_Method(a);
        }

        Dx = Dxc;
        Dy = Dyc;
        Dz = Dzc;

    }

    void Euler_Method(Nucleus[] a) {
        Nucleus[] temp = new Nucleus[a.length];

        double xi;
        double xij, yij, zij, rij, Fijx, Fijy, Fijz;
        int i, j;

        for (i = 0; i < N; i++) {
            temp[i] = a[i];
        }

        for (i = 0; i < N; i++) {

            Fijx = 0.0;
            Fijy = 0.0;
            Fijz = 0.0;

            for (j = 0; j < N; j++) {
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
            temp[i].setX(a[i].getX() - Fijx * DT + Math.sqrt(2.0 * Dx) * Math.sqrt(DT) * xi);

            if (temp[i].getX() < 0) {
                temp[i].setX(0.0);
            }
            if (temp[i].getX() > Lx) {
                temp[i].setX(Lx);
            }

            xi = Box_Muller_Method(0.0, 1.0);
            temp[i].setY(a[i].getY() - Fijy * DT + Math.sqrt(2.0 * Dy) * Math.sqrt(DT) * xi);

            if (temp[i].getY() < 0) {
                temp[i].setY(0.0);
            }
            if (temp[i].getY() > Ly) {
                temp[i].setY(Ly);
            }

            xi = Box_Muller_Method(0.0, 1.0);
            temp[i].setZ(a[i].getZ() - Fijz * DT + Math.sqrt(2.0 * Dz) * Math.sqrt(DT) * xi);

            if (temp[i].getZ() < 0) {
                temp[i].setZ(0.0);
            }
            if (temp[i].getZ() > Lz) {
                temp[i].setZ(Lz);
            }

        }

        for (i = 0; i < N; i++) {
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

    void Tanh_blob(Nucleus[] a, ImageStack output) {
        int nbCPUs = Runtime.getRuntime().availableProcessors();
        IntensThread[] intensThreads = new IntensThread[nbCPUs];
        for (int thread = 0; thread < nbCPUs; thread++) {
            intensThreads[thread] = new IntensThread(thread, nbCPUs, a, output);
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

    void addNoise(ImageStack stack) {
        int nbCPUs = Runtime.getRuntime().availableProcessors();
        NoiseThread[] noiseThreads = new NoiseThread[nbCPUs];
        for (int thread = 0; thread < nbCPUs; thread++) {
            noiseThreads[thread] = new NoiseThread(thread, nbCPUs, stack);
            noiseThreads[thread].start();
        }
        try {
            for (int thread = 0; thread < nbCPUs; thread++) {
                noiseThreads[thread].join();
            }
        } catch (InterruptedException ie) {
            IJ.error("A thread was interrupted during output generation.");
        }
    }

    class IntensThread extends Thread {

        private final Nucleus[] a;
        private final ImageStack output;
        private final int thread;
        private final int nThreads;

        public IntensThread(int thread, int nThreads, Nucleus[] a, ImageStack output) {
            this.thread = thread;
            this.nThreads = nThreads;
            this.a = a;
            this.output = output;
        }

        public void run() {
            for (int k = thread; k < output.getSize(); k += nThreads) {
                for (int j = 0; j < output.getHeight(); j++) {
                    for (int i = 0; i < output.getWidth(); i++) {
                        double xi = px * (2.0 * i + 1) / 2.0;
                        double yi = py * (2.0 * j + 1) / 2.0;
                        double zi = pz * (2.0 * k + 1) / 2.0;
                        double intensity = Iback;
                        double intmax = 0.0;
                        for (int l = 0; l < N; l++) {
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

                            double temp = 1 - (xr * xr / (A * A) + yr * yr / (B * B) + zr * zr / (C * C));
                            temp = (Math.tanh(sigma * temp) + 1.0) / 2.0;
                            temp = (I0 / Math.abs((Math.tanh(sigma) + 1.0) / 2.0)) * temp;

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

    class NoiseThread extends Thread {

        private final ImageStack stack;
        private final int thread;
        private final int nThreads;

        public NoiseThread(int thread, int nThreads, ImageStack stack) {
            this.thread = thread;
            this.nThreads = nThreads;
            this.stack = stack;
        }

        public void run() {
            int width = stack.getWidth();
            int height = stack.getHeight();
            double snr2 = SNR * SNR;
            for (int z = thread + 1; z <= stack.getSize(); z += nThreads) {
                ImageProcessor slice = stack.getProcessor(z);
                for (int y = 0; y < height; y++) {
                    for (int x = 0; x < width; x++) {
                        double I0 = slice.getPixelValue(x, y);
                        double w = I0 / snr2;
                        GammaDistribution dist = new GammaDistribution(snr2, w);
                        slice.putPixelValue(x, y, dist.sample());
                    }
                }
            }
        }
    }
}
