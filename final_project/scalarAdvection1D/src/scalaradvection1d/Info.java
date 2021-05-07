package scalaradvection1d;

/**
 *
 * @author kenny
 */
public class Info {
    public static final int NUM_X_CELLS = 200;
    public static final int NUM_GHOST_CELLS = 3;
    public static final int NUM_RK_STEPS = 3; // can be 1, 2, or 3 here
    public static final int MAX_TIME_ITER = 100000;
    
    public static final double STOPPING_TIME = 2.0;
    public static final double MIN_X = -1.0;
    public static final double MAX_X = 1.0;
    public static final double ADVECTION_VEL = 1.0;
    public static final double CFL = 0.6;
    public static final ReconstructionTypes RECONST_TYPE = ReconstructionTypes.LIMITED_LW;
    
    public static final String BASEPATH = "./output/";
}
