package shallowwater;

/**
 *
 * @author Kangqi Fu
 */
public class Cell {
    public final static int NUM_VARS = 3; // here for this simple equations is three
    
    public double[][] U = new double[Info.NUM_RK_STEPS + 1][NUM_VARS]; // from 1D to 2D to store the system
    public double[] uWest = new double[NUM_VARS];
    public double[] uEast = new double[NUM_VARS];
    public double[] uSouth = new double[NUM_VARS];
    public double[] uNorth = new double[NUM_VARS];
    public double cx;
    public double cy;
    public double dx;
    public double dy;
    public double[] totalFlux = new double[NUM_VARS];
    
}
