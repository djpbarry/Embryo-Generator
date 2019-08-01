/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package image_simulator;

import TimeAndDate.TimeAndDate;
import UtilClasses.GenUtils;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import ij.plugin.Binner;
import ij.plugin.GaussianBlur3D;
import ij.plugin.ImageCalculator;
import ij.plugin.SubstackMaker;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.io.File;
import java.io.IOException;
import java.util.Random;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
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

    private final String TITLE = "Image Simulator";
    double Dx = 0.0, Dy = 0.0, Dz = 0.0; //Diffusion coefficient in x, y, and z directions
    double DT = 0.01; //time grid for SDE

    final double simSizeX = 0.962;
    final double simSizeY = 0.962;
    final double simSizeZ = 0.962; //pixel size in x, y, and z axis

    final double outputSizeX;
    final double outputSizeY;
    final double outputSizeZ;

    double framerate = 1.0; //frame rate of images

    final double Lx = 150;
    final double Ly = 150;
    final double Lz = 150; //domain size

    double I0, sigma = 2.0, Iback = 1000, A = 10.0, B = 7.5, C = 7.5;
    double eps = 1.0, gam = 20.0;
    double distanceThreshold = 90.0;
    int maxframe;

    float PI = 3.14159265358979323846f;

    private Random r = new Random();

    static int count_BM = 0;

    double SNR = 1.5;

    int N = 10 + r.nextInt(25);

    private final String outputDir;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        double px = Double.parseDouble(args[0]);
        double py = Double.parseDouble(args[1]);
        double pz = Double.parseDouble(args[2]);
        (new Image_Simulator(new double[]{px, py, pz}, args[3])).run();
        System.exit(0);
    }

    public Image_Simulator(double[] outputVoxSize, String outputDir) {
        this.outputSizeX = outputVoxSize[0];
        this.outputSizeY = outputVoxSize[1];
        this.outputSizeZ = outputVoxSize[2];
        this.outputDir = outputDir;
    }

    public void run() {
        System.out.println(String.format("%s %s", TITLE, TimeAndDate.getCurrentTimeAndDate()));
        Nucleus[] a = new Nucleus[N];
        int tmax = (int) Math.round(1.5 / DT);

        double Im = 20000; //mean intensity in nucleus

        I0 = Im - Iback;

        //number of pixels in each direction
        int nx = (int) Math.round(Lx / simSizeX);
        int ny = (int) Math.round(Ly / simSizeY);
        int nz = (int) Math.round(Lz / simSizeZ);

        saveNucleiStack(nx, ny, nz, a, tmax);

        saveCellMembraneStack(nx, ny, nz, a);

        IJ.log(String.format("Done %s", TimeAndDate.getCurrentTimeAndDate()));
    }

    void saveNucleiStack(int nx, int ny, int nz, Nucleus[] a, int tmax) {
        ImageStack nucleiOutput = new ImageStack(nx, ny);
        for (int n = 1; n <= nz; n++) {
            nucleiOutput.addSlice(new FloatProcessor(nx, ny));
        }
        maxframe = (int) Math.round(tmax * DT * framerate);

        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Initialising nuclei..."));
        //setting the initial positions of particle
        for (int i = 0; i < N; i++) {
            a[i] = new Nucleus();
            a[i].setTheta_x(0.25 * PI * r.nextFloat());
            a[i].setTheta_y(0.25 * PI * r.nextFloat());
            a[i].setTheta_z(0.5 * PI * r.nextFloat());
            a[i].setX(0.4 * Lx + 0.2 * Lx * r.nextDouble());
            a[i].setY(0.4 * Ly + 0.2 * Ly * r.nextDouble());
            a[i].setZ(0.4 * Lz + 0.2 * Lz * r.nextDouble());

        }
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Simulating nuclei movement..."));
        simulation(a, tmax, nucleiOutput);
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Adding noise to nuclei image..."));
        addNoise(nucleiOutput);
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Saving nuclei image..."));
        try {
            saveStack(new ImagePlus("", nucleiOutput), "Nuclei.tif", nx, ny, nz);
        } catch (Exception e) {
            GenUtils.logError(e, "Encountered problem while saving nuclei images.");
        }
    }

    void saveCellMembraneStack(int nx, int ny, int nz, Nucleus[] a) {
        ImagePlus cellMembraneOutput = generateCellMembraneStack(nx, ny, nz, a);
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Blurring membrane image..."));
        GaussianBlur3D.blur(cellMembraneOutput, 1.0 / simSizeX, 1.0 / simSizeY, 1.0 / simSizeZ);
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Adding noise to membrane image..."));
        addNoise(cellMembraneOutput.getImageStack());
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Saving membrane image..."));
        try {
            saveStack(cellMembraneOutput, "Cell_Membranes.tif", nx, ny, nz);
        } catch (Exception e) {
            GenUtils.logError(e, "Encountered problem while saving cell membrane images.");
        }
    }

    void saveStack(ImagePlus input, String filename, int nx, int ny, int nz) throws DependencyException, ServiceException, FormatException, IOException {
        int size = input.getNSlices();
        int stepZ = (int) Math.round(outputSizeZ / simSizeZ);
        int xBin = (int) Math.round(outputSizeX / simSizeX);
        int yBin = (int) Math.round(outputSizeY / simSizeY);
        ImagePlus subStack = (new SubstackMaker()).makeSubstack(input, String.format("1-%d-%d", size, stepZ));
        input = null;
        ImageStack binnedStack = binStack(subStack.getImageStack(), xBin, yBin);
        ImagePlus imp = new ImagePlus("", binnedStack);
        Calibration cal = imp.getCalibration();
        cal.pixelWidth = simSizeX * nx / binnedStack.getWidth();
        cal.pixelHeight = simSizeY * ny / binnedStack.getHeight();
        cal.pixelDepth = simSizeZ * nz / binnedStack.getSize();
        cal.setXUnit(String.format("%cm", IJ.micronSymbol));
        cal.setYUnit(String.format("%cm", IJ.micronSymbol));
        cal.setZUnit(String.format("%cm", IJ.micronSymbol));
        IJ.saveAs(imp, "TIF", String.format("%s%s%s", outputDir, File.separator, filename));
//        BioFormatsImageWriter.saveStack(binnedStack,
//                new File(String.format("%s%s%s", outputDir, File.separator, filename)),
//                null, FormatTools.FLOAT, "XYZCT", new int[]{binnedStack.getWidth(), binnedStack.getHeight(), binnedStack.getSize(), 1, 1},
//                new double[]{simSizeX * nx / binnedStack.getWidth(),
//                    simSizeY * ny / binnedStack.getHeight(),
//                    simSizeZ * nz / binnedStack.getSize()});
    }

    ImageStack binStack(ImageStack input, int xBin, int yBin) {
        int outWidth = (int) Math.round(input.getWidth() / xBin);
        int outHeight = (int) Math.round(input.getHeight() / yBin);
        ImageStack binnedStack = new ImageStack(outWidth, outHeight);
        Binner b = new Binner();
        for (int s = 1; s <= input.getSize(); s++) {
            ImageProcessor slice = input.getProcessor(s);
            binnedStack.addSlice(b.shrink(slice, xBin, yBin, Binner.SUM));
        }
        return binnedStack;
    }

    ImagePlus generateCellMembraneStack(int nx, int ny, int nz, Nucleus[] a) {
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Generating nuclei centroid image..."));
        ImageHandler pointImage = generatePointImage(nx, ny, nz, a);
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Generating EDT..."));
        ImageHandler edt = EDT.run(pointImage, 1, true, -1);
        ImageHandler edtDup = edt.duplicate();
        float edtMaxValue = (float) edtDup.getMax();
        edtDup.invert();
        edtDup.addValue(edtMaxValue + 1.0f);
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Thresholding EDT..."));
        ImageStack thresholdedEDT = getThresholdedEDT(nx, ny, nz, edtDup, distanceThreshold);
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Outlining EDT..."));
        ImageStack edtOutline = outlineMask(nx, ny, nz, thresholdedEDT);
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Running watershed..."));
        Watershed3D watershed = new Watershed3D(edtDup, pointImage, distanceThreshold, 1);
        ImageInt watershedDams = watershed.getDamImage();
        edtDup = null;
        pointImage = null;
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Summing outline and watershed dams..."));
        ImageStack voronoiPlusMaskOutline = getVoronoiPlusMaskOutline(nx, ny, nz, edtOutline, watershedDams);
        watershedDams = null;
        edtOutline = null;
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Summing watershed dams and EDT..."));
        return getVoronoiPlusEDT(nz, voronoiPlusMaskOutline, edt, thresholdedEDT);
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

    ImageHandler generatePointImage(int nx, int ny, int nz, Nucleus[] a) {
        ImageStack pointOutput = new ImageStack(nx, ny);
        for (int n = 1; n <= nz; n++) {
            ByteProcessor slice = new ByteProcessor(nx, ny);
            slice.setValue(0.0);
            slice.fill();
            pointOutput.addSlice(slice);
        }
        for (int i = 0; i < N; i++) {
            int x = (int) Math.round(a[i].getX() / simSizeX);
            int y = (int) Math.round(a[i].getY() / simSizeY);
            int z = (int) Math.round(a[i].getZ() / simSizeZ);
            pointOutput.setVoxel(x, y, z, 255);
        }
        ImageHandler intImage = ImageInt.wrap(pointOutput);
        intImage.setScale(simSizeX, simSizeZ, "microns");
        return intImage;
    }

    ImageStack getThresholdedEDT(int nx, int ny, int nz, ImageHandler edt, double threshold) {
        ImageStack thresholdedEDT = new ImageStack(nx, ny);
        for (int i = 1; i <= nz; i++) {
            ImageProcessor slice = edt.getImageStack().getProcessor(i);
            slice.setThreshold(0.0, threshold, ImageProcessor.NO_LUT_UPDATE);
            slice = slice.createMask();
            thresholdedEDT.addSlice(slice);
        }
        return thresholdedEDT;
    }

    ImageStack outlineMask(int nx, int ny, int nz, ImageStack input) {
        ImageStack outlinedMask = new ImageStack(nx, ny);
        for (int i = 1; i <= nz; i++) {
            ImageProcessor slice = input.getProcessor(i).duplicate();
            assert (slice instanceof ByteProcessor);
            ((ByteProcessor) slice).outline();
            slice.invert();
            outlinedMask.addSlice(slice);
        }
        return outlinedMask;
    }

    ImageStack getVoronoiPlusMaskOutline(int nx, int ny, int nz, ImageStack thresholdedEDT, ImageHandler lineImage) {
        ImagePlus voronoiPlusMaskOutline = (new ImageCalculator()).run("Add create stack", new ImagePlus("", thresholdedEDT), lineImage.getImagePlus());
        ImageStack voronoiPlusMaskOutlineStack = new ImageStack(nx, ny);
        for (int i = 1; i <= nz; i++) {
            ImageProcessor slice = voronoiPlusMaskOutline.getImageStack().getProcessor(i).convertToFloatProcessor();
            slice.multiply(1.0 / 255.0);
//            slice.add(Iback);
            voronoiPlusMaskOutlineStack.addSlice(slice);
        }
        return voronoiPlusMaskOutlineStack;
    }

    ImagePlus getVoronoiPlusEDT(int nz, ImageStack voronoiPlusMaskOutlineStack, ImageHandler edt, ImageStack mask) {
        ImagePlus voronoiPlusEDT = (new ImageCalculator()).run("Multiply create stack", new ImagePlus("", voronoiPlusMaskOutlineStack), edt.getImagePlus());
        for (int i = 1; i <= nz; i++) {
            ImageProcessor slice = voronoiPlusEDT.getImageStack().getProcessor(i);
            slice.add(1.0);
            slice.multiply(Iback / 3.0);
            ImageProcessor maskSlice = mask.getProcessor(i);
            maskSlice.invert();
            maskSlice.multiply(1.0 / 255.0);
            maskSlice.multiply(1.1);
            maskSlice.add(1.0);
        }
        return (new ImageCalculator()).run("Multiply create stack", voronoiPlusEDT, (new ImagePlus("", mask)));
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
                        double xi = simSizeX * (2.0 * i + 1) / 2.0;
                        double yi = simSizeY * (2.0 * j + 1) / 2.0;
                        double zi = simSizeZ * (2.0 * k + 1) / 2.0;
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
