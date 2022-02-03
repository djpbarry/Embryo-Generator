/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package net.calm.embryogen.image_simulator;

import java.util.Random;

import static java.lang.Math.PI;

/**
 * @author David Barry <david.barry at crick dot ac dot uk>
 */
public class Nucleus {
    //x, y, and z coordinates
    private double x;
    private double y;
    private double z;

    //oriantation around xyz axis
    private double theta_x;
    private double theta_y;
    private double theta_z;

    private double xLength;
    private double yLength;
    private double zLength;

    private final Random r;
    private final double[] initialCoords;
    private final double[] nuclearDimensions;
    private final boolean isCluster;


    public Nucleus(double[] initialCoords, double[] nuclearDimensions, Random r, boolean isCluster) {
        this.r = r;
        this.initialCoords = initialCoords;
        this.nuclearDimensions = nuclearDimensions;
        this.isCluster = isCluster;
        initialise();
    }

    private void initialise() {
        this.setTheta_x(0.25 * PI * r.nextFloat());
        this.setTheta_y(0.25 * PI * r.nextFloat());
        this.setTheta_z(0.5 * PI * r.nextFloat());
        if (isCluster) {
            this.setX(initialCoords[0] * (0.8 + 0.4 * r.nextDouble()));
            this.setY(initialCoords[1] * (0.8 + 0.4 * r.nextDouble()));
            this.setZ(initialCoords[2] * (0.8 + 0.4 * r.nextDouble()));
        } else {
            this.setX(initialCoords[0] + 0.2 * r.nextGaussian());
            this.setY(initialCoords[1] + 0.2 * r.nextGaussian());
            this.setZ(initialCoords[2] + 0.1 * r.nextGaussian());
        }
        this.setxLength(nuclearDimensions[0] + 0.2 * nuclearDimensions[0] * r.nextGaussian());
        this.setyLength(nuclearDimensions[1] + 0.2 * nuclearDimensions[1] * r.nextGaussian());
        this.setzLength(nuclearDimensions[2] + 0.2 * nuclearDimensions[2] * r.nextGaussian());
    }

    public double getX() {
        return x;
    }

    public void setX(double x) {
        this.x = x;
    }

    public double getY() {
        return y;
    }

    public void setY(double y) {
        this.y = y;
    }

    public double getZ() {
        return z;
    }

    public void setZ(double z) {
        this.z = z;
    }

    public double getTheta_x() {
        return theta_x;
    }

    public void setTheta_x(double theta_x) {
        this.theta_x = theta_x;
    }

    public double getTheta_y() {
        return theta_y;
    }

    public void setTheta_y(double theta_y) {
        this.theta_y = theta_y;
    }

    public double getTheta_z() {
        return theta_z;
    }

    public void setTheta_z(double theta_z) {
        this.theta_z = theta_z;
    }

    public double getxLength() {
        return xLength;
    }

    public void setxLength(double xLength) {
        this.xLength = xLength;
    }

    public double getyLength() {
        return yLength;
    }

    public void setyLength(double yLength) {
        this.yLength = yLength;
    }

    public double getzLength() {
        return zLength;
    }

    public void setzLength(double zLength) {
        this.zLength = zLength;
    }

    public double[] getInitialCoords() {
        return initialCoords;
    }

    public double[] getNuclearDimensions() {
        return nuclearDimensions;
    }

    public Random getR() {
        return r;
    }
}
