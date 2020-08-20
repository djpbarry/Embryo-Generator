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
import ij.plugin.filter.Analyzer;
import ij.process.FloatProcessor;
import ij.process.StackStatistics;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.image3d.regionGrowing.Watershed3D;
import net.calm.embryogen.binary.BinaryProcessor;
import net.calm.embryogen.generator.NucleusGenerator;
import net.calm.embryogen.io.StackSaver;
import net.calm.embryogen.noise.NoiseThread;
import net.calm.embryogen.params.SimParams;
import net.calm.embryogen.processor.StackProcessor;
import net.calm.embryogen.simulate.Simulation;
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

    double eps = 1.0, gam = 20.0;
    double distanceThreshold = 75.0;
    double blurRadius = 1.5;
    int maxframe;

    private final int stepZ;
    private final int xBin;
    private final int yBin;

    private final double PI = Math.PI;

    private final Random r = new Random();

    private final double snr;

    //    int N = 10 + r.nextInt(25);
    private final int nCells;

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
        ImageStack nucOutput = generateNucleiStack(nx, ny, nz, a);
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

    ImageStack generateNucleiStack(int nx, int ny, int nz, Nucleus[] a) {
        ImageStack nucleiOutput = new ImageStack(nx, ny);
        for (int n = 1; n <= nz; n++) {
            nucleiOutput.addSlice(new FloatProcessor(nx, ny));
        }
        (new Simulation(params, eps, gam, new double[]{Lx, Ly, Lz})).simulation(a, nucleiOutput);
        return nucleiOutput;
    }

    ImagePlus generateCellMembraneStack(int nx, int ny, int nz, Nucleus[] a) {
        BinaryProcessor bp = new BinaryProcessor(params);
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Generating nuclei centroid image..."));
        ImageHandler pointImage = bp.generatePointImage(nx, ny, nz, a);
//        saveStack(pointImage.getImagePlus(), "point_image.tif");
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Generating EDT..."));
        ImageHandler edt = EDT.run(pointImage, 1, true, -1);
//        saveStack(edt.getImagePlus(), "edt_image.tif");
        ImageHandler edtDup = edt.duplicate();
        float edtMaxValue = (float) edtDup.getMax();
        edtDup.invert();
        edtDup.addValue(edtMaxValue + 1.0f);
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Thresholding EDT..."));
        ImageStack thresholdedEDT = bp.getThresholdedEDT(nx, ny, nz, edtDup, distanceThreshold);
//        saveStack(new ImagePlus("", thresholdedEDT), "thresh_edt_image.tif");
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Outlining EDT..."));
        ImageStack edtOutline = bp.outlineMask(nx, ny, nz, thresholdedEDT);
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
        ImageStack voronoiPlusMaskOutline = bp.getVoronoiPlusMaskOutline(nx, ny, nz, edtOutline, watershedDams);
//        saveStack(new ImagePlus("", voronoiPlusMaskOutline), "dams_plus_outline_image.tif");
        watershedDams = null;
        edtOutline = null;
        System.out.println(String.format("%s %s", TimeAndDate.getCurrentTimeAndDate(), "Summing watershed dams and EDT..."));
        return bp.getVoronoiPlusEDT(nz, voronoiPlusMaskOutline, edt, thresholdedEDT);
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
