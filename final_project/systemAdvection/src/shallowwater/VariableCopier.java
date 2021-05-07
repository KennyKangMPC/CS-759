package shallowwater;

/**
 *
 * @author Kangqi Fu
 */
public class VariableCopier {
    private VariableCopier(){
    }
    
    public static void copyToZerothRKStep(Cell[][] cells){
        for (int i = Info.NUM_GHOST_CELLS; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS; i++){
            for (int j = Info.NUM_GHOST_CELLS; j < Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS; j++){
                System.arraycopy(cells[i][j].U[Info.NUM_RK_STEPS], 0, cells[i][j].U[0], 0, Cell.NUM_VARS); //This might have error check back :)
            }
        }
    }
}
