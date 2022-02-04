package net.calm.embryogen.image_simulator;

import java.util.ArrayList;
import java.util.Random;

public class NucleusGroup {
    final ArrayList<Nucleus> children;
    final Nucleus parent;

    public NucleusGroup(Nucleus parent, int nChildren) {
        this.parent = parent;
        this.children = new ArrayList<>();
        initiateChildren(nChildren);
    }

    private void initiateChildren(int nChildren) {
        double[] parentNucDims = parent.getNuclearDimensions();
        Random r = parent.getR();
        double[] childNucDims = new double[]{
                parentNucDims[0] * 0.2 + 0.1 * r.nextDouble(),
                parentNucDims[1] * 0.2 + 0.1 * r.nextDouble(),
                parentNucDims[2] * 0.2 + 0.1 * r.nextDouble()
        };
        for (int i = 0; i < nChildren; i++) {
            double[] childInitCoords = new double[]{
                    parent.getX() + Math.max(parentNucDims[0], parentNucDims[1]) * r.nextGaussian(),
                    parent.getY() + Math.max(parentNucDims[0], parentNucDims[1]) * r.nextGaussian(),
                    parent.getZ() + 0.5 * parentNucDims[2] * r.nextGaussian()};
            children.add(
                    new Nucleus(childInitCoords, childNucDims, parent.getR(), Nucleus.FIXED));
        }
    }

    public ArrayList<Nucleus> getChildren() {
        return children;
    }

    public ArrayList<Nucleus> getAllElements() {
        ArrayList<Nucleus> allElements = new ArrayList<>();
        allElements.add(getParent());
        allElements.addAll(getChildren());
        return allElements;
    }

    public Nucleus getParent() {
        return parent;
    }
}
