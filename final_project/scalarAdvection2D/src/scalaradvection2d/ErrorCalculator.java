package scalaradvection2d;

/**
 *
 * @author Kangqi Fu
 */
public class ErrorCalculator { // This can be easily transformed into a 2D form
    
    double errorL1, errorL2, errorLinf;
    
    public ErrorCalculator(Cell[][] cells){
        //get the initial distribution
        Cell[][] initCells = new Cell[cells.length][cells[0].length];
        
        for (int i = 0; i < initCells.length; i++){
            for (int j = 0; j < initCells[0].length; j++){
                initCells[i][j] = new Cell();
                initCells[i][j].cx = cells[i][j].cx;
                initCells[i][j].cy = cells[i][j].cy;
                initCells[i][j].dx = cells[i][j].dx;
                initCells[i][j].dy = cells[i][j].dy;
            }

        }
        SchemeInitializer.initializeSolution(initCells);
        
        errorL1 = 0.0;
        errorL2 = 0.0;
        errorLinf = 0.0;
        
        for (int i = Info.NUM_GHOST_CELLS; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS; i++){
            for (int j = Info.NUM_GHOST_CELLS; j < Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS; j++) {
                double error = Math.abs(initCells[i][j].u[0] - cells[i][j].u[0]);
                errorL1 += error;
                errorL2 += error * error;
                errorLinf = Math.max(error, errorLinf);
            }
        }
        int numCells = Info.NUM_X_CELLS * Info.NUM_Y_CELLS;
        errorL1 /= numCells;
        errorL2 = Math.sqrt(errorL2 / numCells);
    }
    
    @Override
    public String toString(){
        return String.format("%-25s%-25s%-25s\n" + "%-25.15e%-25.15e%-25.15e\n", "L1_Error", "L2_Error", "Linf_Error", errorL1, errorL2, errorLinf);
    }
    
}
