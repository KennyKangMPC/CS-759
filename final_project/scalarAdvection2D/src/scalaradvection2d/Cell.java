package scalaradvection2d;

/**
 *
 * @author Kangqi Fu
 */
public class Cell {
    public double[] u = new double[Info.NUM_RK_STEPS + 1];
    public double uWest, uEast, uSouth, uNorth; // crosponding to u_(i-1/2) and u_(i+1/2)
    public double cx;
    public double cy;
    public double dx;
    public double dy;
    public double totalFlux;
}
