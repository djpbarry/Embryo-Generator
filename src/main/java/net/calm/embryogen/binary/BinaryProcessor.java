package net.calm.embryogen.binary;

import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.ImageCalculator;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import net.calm.embryogen.image_simulator.Nucleus;
import net.calm.embryogen.image_simulator.NucleusGroup;
import net.calm.embryogen.params.SimParams;

public class BinaryProcessor {
    private final SimParams params;

    public BinaryProcessor(SimParams params) {
        this.params = params;
    }

    public ImageHandler generatePointImage(int nx, int ny, int nz, NucleusGroup[] a) {
        ImageStack pointOutput = new ImageStack(nx, ny);
        for (int n = 1; n <= nz; n++) {
            ByteProcessor slice = new ByteProcessor(nx, ny);
            slice.setValue(0.0);
            slice.fill();
            pointOutput.addSlice(slice);
        }
        for (int i = 0; i < a.length; i++) {
            int x = (int) Math.round(a[i].getParent().getX() / params.getSimSizeX());
            int y = (int) Math.round(a[i].getParent().getY() / params.getSimSizeY());
            int z = (int) Math.round(a[i].getParent().getZ() / params.getSimSizeZ());
            pointOutput.setVoxel(x, y, z, 255);
        }
        ImageHandler intImage = ImageInt.wrap(pointOutput);
        intImage.setScale(params.getSimSizeX(), params.getSimSizeZ(), "microns");
        return intImage;
    }

    public ImageStack getThresholdedEDT(int nx, int ny, int nz, ImageHandler edt, double threshold) {
        ImageStack thresholdedEDT = new ImageStack(nx, ny);
        for (int i = 1; i <= nz; i++) {
            ImageProcessor slice = edt.getImageStack().getProcessor(i);
            slice.setThreshold(0.0, threshold, ImageProcessor.NO_LUT_UPDATE);
            slice = slice.createMask();
            thresholdedEDT.addSlice(slice);
        }
        return thresholdedEDT;
    }

    public ImageStack outlineMask(int nx, int ny, int nz, ImageStack input) {
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

    public ImageStack getVoronoiPlusMaskOutline(int nx, int ny, int nz, ImageStack thresholdedEDT, ImageHandler lineImage) {
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

    public ImagePlus getVoronoiPlusEDT(int nz, ImageStack voronoiPlusMaskOutlineStack, ImageHandler edt, ImageStack mask) {
        ImagePlus voronoiPlusEDT = (new ImageCalculator()).run("Multiply create stack", new ImagePlus("", voronoiPlusMaskOutlineStack), edt.getImagePlus());
        ImageStack scaledMask = new ImageStack(mask.getWidth(), mask.getHeight());
        for (int i = 1; i <= nz; i++) {
            double intensScale = params.getMembraneMaxIntensScale() - (params.getMembraneMaxIntensScale() - params.getMembraneMinIntensScale()) * (i - 1) / nz;
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
}
