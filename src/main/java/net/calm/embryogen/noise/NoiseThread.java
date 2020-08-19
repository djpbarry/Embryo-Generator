package net.calm.embryogen.noise;

import ij.ImageStack;
import ij.process.ImageProcessor;
import org.apache.commons.math3.distribution.GammaDistribution;

public class NoiseThread extends Thread {

    private final ImageStack stack;
    private final int thread;
    private final int nThreads;
    private final double snr;

    public NoiseThread(int thread, int nThreads, double snr, ImageStack stack) {
        this.thread = thread;
        this.nThreads = nThreads;
        this.stack = stack;
        this.snr = snr;
    }

    public void run() {
        int width = stack.getWidth();
        int height = stack.getHeight();
        double snr2 = Math.pow(snr, 2.0);
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
