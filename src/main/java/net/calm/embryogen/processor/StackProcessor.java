package net.calm.embryogen.processor;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import ij.plugin.Binner;
import ij.plugin.SubstackMaker;
import ij.process.ImageProcessor;
import net.calm.embryogen.params.SimParams;

public class StackProcessor {

    public static ImagePlus downsizeStack(ImagePlus input, int nx, int ny, int nz, int xBin, int yBin, int stepZ, SimParams params) {
        ImagePlus subStack = downSizeStack(input, stepZ);
        input = null;
        ImageStack binnedStack = binStack(subStack.getImageStack(), xBin, yBin, Binner.SUM);
        ImagePlus imp = new ImagePlus("", binnedStack);
        Calibration cal = imp.getCalibration();
        cal.pixelWidth = params.getSimSizeX() * nx / binnedStack.getWidth();
        cal.pixelHeight = params.getSimSizeY() * ny / binnedStack.getHeight();
        cal.pixelDepth = params.getSimSizeZ() * nz / binnedStack.getSize();
        cal.setXUnit(String.format("%cm", IJ.micronSymbol));
        cal.setYUnit(String.format("%cm", IJ.micronSymbol));
        cal.setZUnit(String.format("%cm", IJ.micronSymbol));
        return imp;
    }

    public static ImagePlus downSizeStack(ImagePlus input, int stepZ) {
        int size = input.getNSlices();
        return (new SubstackMaker()).makeSubstack(input, String.format("1-%d-%d", size, stepZ));
    }

    public static ImageStack binStack(ImageStack input, int xBin, int yBin, int binningMethod) {
        int outWidth = (int) Math.round(input.getWidth() / xBin);
        int outHeight = (int) Math.round(input.getHeight() / yBin);
        ImageStack binnedStack = new ImageStack(outWidth, outHeight);
        Binner b = new Binner();
        for (int s = 1; s <= input.getSize(); s++) {
            ImageProcessor slice = input.getProcessor(s);
            binnedStack.addSlice(b.shrink(slice, xBin, yBin, binningMethod));
        }
        return binnedStack;
    }
}

