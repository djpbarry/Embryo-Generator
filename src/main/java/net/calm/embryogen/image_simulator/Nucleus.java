/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package image_simulator;

/**
 *
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

}
