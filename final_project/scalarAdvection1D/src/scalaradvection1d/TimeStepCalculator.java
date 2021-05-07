package scalaradvection1d;

/**
 *
 * @author kenny
 */
public class TimeStepCalculator {
    private TimeStepCalculator(){
    }
    
    public static double getTimeStep(Cell[] cells){
        // Assumed that dx is constant for all cells
        return cells[0].dx / Math.abs(Info.ADVECTION_VEL) * Info.CFL; 
    }
}
