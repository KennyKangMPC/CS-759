package shallowwater;

/**
 *
 * @author kenny
 */
public class TimeStepCalculator {
    //TODO: need to generalize this function make it more general
    private TimeStepCalculator(){
    }
    
    public static double getTimeStep(Cell[][] cells){
        
        double minDT = Double.POSITIVE_INFINITY;
        
        for (int i = Info.NUM_GHOST_CELLS; i < cells.length - Info.NUM_GHOST_CELLS; i++){
            for (int j = Info.NUM_GHOST_CELLS; j < cells[i].length - Info.NUM_GHOST_CELLS; j++){
                double cell_dt = getTimeStepForCell(cells[i][j]);
                if (minDT > cell_dt)
                    minDT = cell_dt;
            }
        }
        return minDT * Info.CFL;
    }
    
    private static double getTimeStepForCell(Cell cell){
        int rkStep = 0;
        double h = cell.U[rkStep][0];
        double u = cell.U[rkStep][1] / h;
        double v = cell.U[rkStep][2] / h;
        double gravity = 9.81;
        double c = Math.sqrt(gravity * h);
        double dtx = cell.dx / (Math.abs(u) + c);
        double dty = cell.dy / (Math.abs(v) + c);
        return 1.0 / (1.0 / dtx + 1.0 / dty);
    }
}
