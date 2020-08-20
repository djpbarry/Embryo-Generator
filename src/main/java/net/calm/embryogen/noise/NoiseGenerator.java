package net.calm.embryogen.noise;

import ij.IJ;
import ij.ImageStack;

public class NoiseGenerator {
    public static void addNoise(ImageStack stack, double snr) {
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
