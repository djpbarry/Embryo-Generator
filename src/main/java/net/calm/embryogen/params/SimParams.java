package net.calm.embryogen.params;

public class SimParams {
    private double simSizeX;
    private double simSizeY;
    private double simSizeZ; //pixel size in x, y, and z axis
    private double outputSizeX;
    private double outputSizeY;
    private double outputSizeZ;
    private final double sigma = 2.0, Iback = 100, A = 10.0, B = 7.5, C = 7.5;
    private final double nucMaxIntens = 1000.0; //mean intensity in nucleus
    private final double nucMinIntens = 0.0;
    private double Dx = 0.0, Dy = 0.0, Dz = 0.0; //Diffusion coefficient in x, y, and z directions
    private final double DT = 0.01; //time grid for SDE
    private final double framerate = 1.0; //frame rate of images
    private final double membraneMinIntensScale = 0.0;
    private final double membraneMaxIntensScale = 1.0;

    public SimParams() {

    }

    public double getSimSizeX() {
        return simSizeX;
    }

    public void setSimSizeX(double simSizeX) {
        this.simSizeX = simSizeX;
    }

    public double getSimSizeY() {
        return simSizeY;
    }

    public void setSimSizeY(double simSizeY) {
        this.simSizeY = simSizeY;
    }

    public double getSimSizeZ() {
        return simSizeZ;
    }

    public void setSimSizeZ(double simSizeZ) {
        this.simSizeZ = simSizeZ;
    }

    public double getOutputSizeX() {
        return outputSizeX;
    }

    public void setOutputSizeX(double outputSizeX) {
        this.outputSizeX = outputSizeX;
    }

    public double getOutputSizeY() {
        return outputSizeY;
    }

    public void setOutputSizeY(double outputSizeY) {
        this.outputSizeY = outputSizeY;
    }

    public double getOutputSizeZ() {
        return outputSizeZ;
    }

    public void setOutputSizeZ(double outputSizeZ) {
        this.outputSizeZ = outputSizeZ;
    }

    public double getSigma() {
        return sigma;
    }

    public double getIback() {
        return Iback;
    }

    public double getA() {
        return A;
    }

    public double getB() {
        return B;
    }

    public double getC() {
        return C;
    }

    public double getNucMaxIntens() {
        return nucMaxIntens;
    }

    public double getNucMinIntens() {
        return nucMinIntens;
    }

    public double getDx() {
        return Dx;
    }

    public double getDy() {
        return Dy;
    }

    public double getDz() {
        return Dz;
    }

    public double getDT() {
        return DT;
    }

    public double getFramerate() {
        return framerate;
    }

    public void setDx(double dx) {
        Dx = dx;
    }

    public void setDy(double dy) {
        Dy = dy;
    }

    public void setDz(double dz) {
        Dz = dz;
    }

    public double getMembraneMinIntensScale() {
        return membraneMinIntensScale;
    }

    public double getMembraneMaxIntensScale() {
        return membraneMaxIntensScale;
    }
}
