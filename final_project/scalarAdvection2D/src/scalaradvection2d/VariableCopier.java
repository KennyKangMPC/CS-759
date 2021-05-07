package scalaradvection2d;

/**
 *
 * @author Kangqi Fu
 */
public class VariableCopier {
    private VariableCopier(){
    }
    
    public static void copyToZerothRKStep(Cell[][] cells){
        for (int i = Info.NUM_GHOST_CELLS; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS; i++){
            for (int j = Info.NUM_GHOST_CELLS; j < Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS; j++)
                cells[i][j].u[0] = cells[i][j].u[Info.NUM_RK_STEPS];
        }
    }
}
