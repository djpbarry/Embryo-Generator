/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package net.calm.embryogen.image_simulator;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.ResultsTable;
import ij.plugin.Binner;
import ij.plugin.GaussianBlur3D;
import ij.plugin.ImageCalculator;
import ij.plugin.filter.Analyzer;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.StackStatistics;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.image3d.regionGrowing.Watershed3D;
import net.calm.embryogen.generator.NucleusGenerator;
import net.calm.embryogen.intensity.IntensThread;
import net.calm.embryogen.io.StackSaver;
import net.calm.embryogen.noise.NoiseThread;
import net.calm.embryogen.params.SimParams;
import net.calm.embryogen.processor.StackProcessor;
import net.calm.iaclasslibrary.IO.DataWriter;
import net.calm.iaclasslibrary.TimeAndDate.TimeAndDate;
import net.calm.iaclasslibrary.UtilClasses.GenUtils;

import java.io.File;
import java.io.IOException;
import java.util.Random;

/**
 * @author David Barry <david.barry at crick dot ac dot uk>
 */
public class Image_Simulator {

    private final String TITLE = "Image Simulator";

    final double Lx = 150;
    final double Ly = 150;
    final double Lz = 0.5; //domain size

    double membraneMinIntensScale = 0.0;
    double membraneMaxIntensScale = 1.0;
    double eps = 1.0, gam = 20.0;
    double distanceThreshold = 75.0;
    double blurRadius = 1.5;
    int maxframe;

    private final int stepZ;
    private final int xBin;
    private final int yBin;

    private final double PI = Math.PI;

    private Random r = new Random();

    static int count_BM = 0;

    private final double snr;

    //    int N = 10 + r.nextInt(25);
    private int nCells;

    private final String simOutputDir;

    private final String groundTruthDir;

    private final ResultsTable resultsTable = Analyzer.getResultsTable();

    private final SimParams params;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        double px = Double.parseDouble(args[0]);
        double py = Double.parseDouble(args[1]);
        double pz = Double.parseDouble(args[2]);
        double sx = Double.parseDouble(args[3]);
        double sy = Double.parseDouble(args[4]);
        double sz = Double.parseDouble(args[5]);
        System.setProperty("java.awt.headless", "true");
        (new Image_Simulator(new double[]{px, py, pz}, new double[]{sx, sy, sz}, Double.parseDouble(args[6]), Integer.parseInt(args[7]), args[8], args[9])).run();
        System.exit(0);
    }

    public Image_Simulator(double[] outputVoxSize, double[] simVoxSize, double snr, int nCells, String simOutputDir, String groundTruthOutputDir) {
        this.params = new SimParams();
        params.setOutputSizeX(outputVoxSize[0]);
        params.setOutputSizeY(outputVoxSize[1]);
        params.setOutputSizeZ(outputVoxSize[2]);
        params.setSimSizeX(simVoxSize[0]);
        params.setSimSizeY(simVoxSize[1]);
        params.setSimSizeZ(simVoxSize[2]);
        this.snr = snr;
        this.nCells = nCells;
        this.simOutputDir = GenUtils.openResultsDirectory(String.format("%s%s%s_snr%f_ncells%d", simOutputDir, File.separator, TITLE, this.snr, this.nCells));
        this.groundTruthDir = GenUtils.openResultsDirectory(String.format("%s%s%s_snr%f_ncells%d", groundTruthOutputDir, File.separator, TITLE, this.snr, this.nCells));
        this.stepZ = (int) Math.round(params.getOutputSizeX() / params.getSimSizeX());
        this.xBin = (int) Math.round(params.getOutputSizeY() / params.getSimSizeY());
        this.yBin = (int) Math.round(params.getOutputSizeZ() / params.getSimSizeZ());
    }

    public void run() {
        System.out.println(String.format("%s %s", TITLE, TimeAndDate.getCurrentTimeAndDate()));
        System.out.println(String.format("Simulation_Size_X = %f", params.getSimSizeX()));
        System.out.println(String.format("Simulation_Size_Y = %f", params.getSimSizeY()));
        System.out.println(String.format("Simulation_Size_Z = %f", params.getSimSizeZ()));
        System.out.println(String.format("Output_Size_X = %f", params.getOutputSizeX()));
        System.out.println(String.format("Output_Size_Y = %f", params.getOutputSizeY()));
        System.out.println(String.format("Output_Size_Z = %f", params.getOutputSizeZ()));
        System.out.println(String.format("SNR = %f", snr));
        Nucleus[] a = new Nucleus[nCells];
        int tmax = (int) Math.round(1.5 / params.getDT());
        maxframe = (int) Math.round(tmax * params.getDT() * params.getFramerate());

        //number of pixels in each direction
        int nx = (int) Math.round(Lx / params.getSimSizeX());
        int ny = (int) Math.round(Ly / params.getSimSizeY());
        int nz = (int) Math.round(Lz / params.getSimSizeZ());

        (new NucleusGenerator(nCells, new double[]{Lx, Ly, Lz}, params)).initialiseNuclei(a);
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Simulating nuclei movement..."));
        ImageStack nucOutput = generateNucleiStack(nx, ny, nz, a, tmax);
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Downsizing nuclei image..."));
        ImagePlus downSizedNucleiImage = StackProcessor.downsizeStack(new ImagePlus("", nucOutput), nx, ny, nz, xBin, yBin, stepZ, params);
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Adding noise to nuclei image..."));
        addNoise(downSizedNucleiImage.getImageStack());
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Saving nuclei image..."));
        StackSaver.saveNucleiStack(downSizedNucleiImage, simOutputDir);

        for (int i = 0; i < a.length; i++) {
            resultsTable.setValue("Nucleus_Centroid_X", i, a[i].getX());
            resultsTable.setValue("Nucleus_Centroid_Y", i, a[i].getY());
            resultsTable.setValue("Nucleus_Centroid_Z", i, a[i].getZ());
        }

        ImagePlus cellMembraneOutput = generateCellMembraneStack(nx, ny, nz, a);
//        saveStack(cellMembraneOutput, "raw_cell_membrane_image.tif");
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Blurring membrane image..."));
        GaussianBlur3D.blur(cellMembraneOutput, blurRadius / params.getSimSizeX(), blurRadius / params.getSimSizeY(), blurRadius / params.getSimSizeZ());
//        saveStack(cellMembraneOutput, "blurred_cell_membrane_image.tif");
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Downsizing membrane image..."));
        ImagePlus downSizedCellMembraneImage = StackProcessor.downsizeStack(cellMembraneOutput, nx, ny, nz, xBin, yBin, stepZ, params);
//        saveStack(downSizedCellMembraneImage, "downsized_cell_membrane_image.tif");
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Adding noise to membrane image..."));
        addNoise(downSizedCellMembraneImage.getImageStack());
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Saving membrane image..."));
        StackSaver.saveCompositeStack(downSizedCellMembraneImage, simOutputDir, String.format("Sim_Image_snr%f_ncells%d.tif", snr, nCells));

        try {
            DataWriter.saveResultsTable(resultsTable, new File(String.format("%s%s%s_snr%f_ncells%d.csv", groundTruthDir, File.separator, "Ground_Truth_Data", snr, nCells)));
        } catch (IOException e) {
            GenUtils.logError(e, "Error saving ground truth data.");
        }

        IJ.log(String.format("Done %s", TimeAndDate.getCurrentTimeAndDate()));
    }

    ImageStack generateNucleiStack(int nx, int ny, int nz, Nucleus[] a, int tmax) {
        ImageStack nucleiOutput = new ImageStack(nx, ny);
        for (int n = 1; n <= nz; n++) {
            nucleiOutput.addSlice(new FloatProcessor(nx, ny));
        }
        simulation(a, tmax, nucleiOutput);
        return nucleiOutput;
    }

    ImagePlus generateCellMembraneStack(int nx, int ny, int nz, Nucleus[] a) {
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Generating nuclei centroid image..."));
        ImageHandler pointImage = generatePointImage(nx, ny, nz, a);
//        saveStack(pointImage.getImagePlus(), "point_image.tif");
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Generating EDT..."));
        ImageHandler edt = EDT.run(pointImage, 1, true, -1);
//        saveStack(edt.getImagePlus(), "edt_image.tif");
        ImageHandler edtDup = edt.duplicate();
        float edtMaxValue = (float) edtDup.getMax();
        edtDup.invert();
        edtDup.addValue(edtMaxValue + 1.0f);
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Thresholding EDT..."));
        ImageStack thresholdedEDT = getThresholdedEDT(nx, ny, nz, edtDup, distanceThreshold);
//        saveStack(new ImagePlus("", thresholdedEDT), "thresh_edt_image.tif");
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Outlining EDT..."));
        ImageStack edtOutline = outlineMask(nx, ny, nz, thresholdedEDT);
//        saveStack(new ImagePlus("", edtOutline), "outlined_edt_image.tif");
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Running watershed..."));
        Watershed3D watershed = new Watershed3D(edtDup, pointImage, distanceThreshold, 1);
        ImageHandler watershedDams = watershed.getDamImage();
        ImagePlus subStack = StackProcessor.downSizeStack(watershed.getWatershedImage3D().getImagePlus(), stepZ);
        ImageStack binnedStack = StackProcessor.binStack(subStack.getImageStack(), xBin, yBin, Binner.MIN);
        ImagePlus binnedImp = new ImagePlus("", binnedStack);
        printGroundTruthResults(binnedStack, nx, ny, nz, a);
        StackSaver.saveGroundTruth(binnedImp, groundTruthDir);
        edtDup = null;
        pointImage = null;
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Summing outline and watershed dams..."));
        ImageStack voronoiPlusMaskOutline = getVoronoiPlusMaskOutline(nx, ny, nz, edtOutline, watershedDams);
//        saveStack(new ImagePlus("", voronoiPlusMaskOutline), "dams_plus_outline_image.tif");
        watershedDams = null;
        edtOutline = null;
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Summing watershed dams and EDT..."));
        return getVoronoiPlusEDT(nz, voronoiPlusMaskOutline, edt, thresholdedEDT);
    }

    void printGroundTruthResults(ImageStack binnedStack, int nx, int ny, int nz, Nucleus[] a) {
        StackStatistics stats = new StackStatistics(new ImagePlus("", binnedStack));
        int[] hist = stats.histogram16;
        double binnedXRes = params.getSimSizeX() * nx / binnedStack.getWidth();
        double binnedYRes = params.getSimSizeY() * ny / binnedStack.getHeight();
        double binnedZRes = params.getSimSizeZ() * nz / binnedStack.getSize();
        for (int i = 0; i < a.length; i++) {
            int x = (int) Math.round(a[i].getX() / binnedXRes);
            int y = (int) Math.round(a[i].getY() / binnedYRes);
            int z = (int) Math.round(a[i].getZ() / binnedZRes);
            if (z >= binnedStack.getSize()) z = binnedStack.getSize() - 1;
            int index = (int) Math.round(binnedStack.getVoxel(x, y, z));
            double vol = hist[index] * binnedXRes * binnedYRes * binnedZRes;
            resultsTable.setValue("Cell_Index", i, index - 1);
            resultsTable.setValue("Cell_Volume_Microns_Cubed", i, vol);
        }
    }

    void simulation(Nucleus[] a, int tmax, ImageStack output) {
        initial_relaxation(a);
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

    void Euler_Method(Nucleus[] a) {
        Nucleus[] temp = new Nucleus[a.length];

        double xi;
        double xij, yij, zij, rij, Fijx, Fijy, Fijz;
        int i, j;

        for (i = 0; i < nCells; i++) {
            temp[i] = a[i];
        }

        for (i = 0; i < nCells; i++) {

            Fijx = 0.0;
            Fijy = 0.0;
            Fijz = 0.0;

            for (j = 0; j < nCells; j++) {
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

        for (i = 0; i < nCells; i++) {
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

    ImageHandler generatePointImage(int nx, int ny, int nz, Nucleus[] a) {
        ImageStack pointOutput = new ImageStack(nx, ny);
        for (int n = 1; n <= nz; n++) {
            ByteProcessor slice = new ByteProcessor(nx, ny);
            slice.setValue(0.0);
            slice.fill();
            pointOutput.addSlice(slice);
        }
        for (int i = 0; i < nCells; i++) {
            int x = (int) Math.round(a[i].getX() / params.getSimSizeX());
            int y = (int) Math.round(a[i].getY() / params.getSimSizeY());
            int z = (int) Math.round(a[i].getZ() / params.getSimSizeZ());
            pointOutput.setVoxel(x, y, z, 255);
        }
        ImageHandler intImage = ImageInt.wrap(pointOutput);
        intImage.setScale(params.getSimSizeX(), params.getSimSizeZ(), "microns");
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
        ImageStack scaledMask = new ImageStack(mask.getWidth(), mask.getHeight());
        for (int i = 1; i <= nz; i++) {
            double intensScale = membraneMaxIntensScale - (membraneMaxIntensScale - membraneMinIntensScale) * (i - 1) / nz;
//            System.out.println(String.format("z: %d IntensScale: %f",i, intensScale));
            ImageProcessor slice = voronoiPlusEDT.getImageStack().getProcessor(i);
            slice.add(1.0);
            slice.multiply(params.getIback());
            ImageProcessor maskSlice = mask.getProcessor(i);
            maskSlice.invert();
            ImageProcessor floatMaskSlice = maskSlice.convertToFloatProcessor();
            floatMaskSlice.multiply(intensScale / 255.0);
            floatMaskSlice.add(1.0);
            scaledMask.addSlice(floatMaskSlice);
        }
        return (new ImageCalculator()).run("Multiply create stack", voronoiPlusEDT, (new ImagePlus("", scaledMask)));
    }

    void addNoise(ImageStack stack) {
        int nbCPUs = Runtime.getRuntime().availableProcessors();
        NoiseThread[] noiseThreads = new NoiseThread[nbCPUs];
        for (int thread = 0; thread < nbCPUs; thread++) {
            noiseThreads[thread] = new NoiseThread(thread, nbCPUs, snr, stack);
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

}
