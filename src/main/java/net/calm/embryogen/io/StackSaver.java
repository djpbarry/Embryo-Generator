package net.calm.embryogen.io;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.Concatenator;
import ij.plugin.HyperStackConverter;
import net.calm.iaclasslibrary.UtilClasses.GenUtils;

import java.io.File;

public class StackSaver {
    public static void saveNucleiStack(ImagePlus nucleiOutput, String simOutputDir) {
        try {
            saveStack(nucleiOutput, "Nuclei.tif", simOutputDir);
        } catch (Exception e) {
            GenUtils.logError(e, "Encountered problem while saving nuclei images.");
        }
    }

    public static void saveStack(ImagePlus input, String filename, String simOutputDir) {
        IJ.saveAs(input, "TIF", String.format("%s%s%s", simOutputDir, File.separator, filename));
//        BioFormatsImageWriter.saveStack(binnedStack,
//                new File(String.format("%s%s%s", outputDir, File.separator, filename)),
//                null, FormatTools.FLOAT, "XYZCT", new int[]{binnedStack.getWidth(), binnedStack.getHeight(), binnedStack.getSize(), 1, 1},
//                new double[]{params.getSimSizeX() * nx / binnedStack.getWidth(),
//                    params.getSimSizeY() * ny / binnedStack.getHeight(),
//                    params.getSimSizeZ() * nz / binnedStack.getSize()});
    }

    public static void saveCompositeStack(ImagePlus cellMembraneOutput, String simOutputDir, String imageFileName) {
        try {
            String nucFileName = String.format("%s%s%s", simOutputDir, File.separator, "Nuclei.tif");
            ImagePlus nuclei = IJ.openImage(nucFileName);
            ImagePlus concat = (new Concatenator()).concatenate(new ImagePlus[]{nuclei, cellMembraneOutput}, false);
            StackSaver.saveStack(HyperStackConverter.toHyperStack(concat, 2, concat.getNSlices(), 1, "xyzct", "composite"), imageFileName, simOutputDir);
            nuclei.changes = false;
            nuclei.close();
            (new File(String.format("%s%s%s", simOutputDir, File.separator, "Nuclei.tif"))).delete();
        } catch (Exception e) {
            GenUtils.logError(e, "Encountered problem while saving cell membrane images.");
        }
    }
}
