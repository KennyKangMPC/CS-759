package scalaradvection2d;

/**
 *
 * @author kenny
 */
public class TimeStepCalculator {
    private TimeStepCalculator(){
    }
    
    public static double getTimeStep(Cell[][] cells){
        // Assumed that dx is constant for all cells
        double dtx = cells[0][0].dx / Math.abs(Info.ADVECTION_VEL_X);
        double dty = cells[0][0].dy / Math.abs(Info.ADVECTION_VEL_Y);
        double dt = 1.0 / (1.0 / dtx + 1.0 / dty);
        return dt * Info.CFL; 
    }
}
