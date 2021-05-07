package scalaradvection1d;

/**
 * Function of 
 * @author Kangqi Fu
 */
public class ErrorCalculator { // This can be easily transformed into a 2D form
    
    double errorL1, errorL2, errorLinf;
    
    public ErrorCalculator(Cell[] cells){
        //get the initial distribution
        Cell[] initCells = new Cell[cells.length];
        
        for (int i = 0; i < initCells.length; i++){
            initCells[i] = new Cell();
            initCells[i].cx = cells[i].cx;
            initCells[i].dx = cells[i].dx;
        }
        SchemeInitializer.initializeSolution(initCells);
        
        errorL1 = 0.0;
        errorL2 = 0.0;
        errorLinf = 0.0;
        
        for (int i = Info.NUM_GHOST_CELLS; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS; i++){
            double error = Math.abs(initCells[i].u[0] - cells[i].u[0]);
            errorL1 += error;
            errorL2 += error * error;
            errorLinf = Math.max(error, errorLinf);
        }
        errorL1 /= Info.NUM_X_CELLS;
        errorL2 = Math.sqrt(errorL2 / Info.NUM_X_CELLS);
    }
    
    @Override
    public String toString(){
        return String.format("%-25s%-25s%-25s\n" + "%-25.15e%-25.15e%-25.15e\n", "L1_Error", "L2_Error", "Linf_Error", errorL1, errorL2, errorLinf);
    }
    
}
