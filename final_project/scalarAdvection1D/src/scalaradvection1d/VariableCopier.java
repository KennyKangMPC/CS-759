package scalaradvection1d;

/**
 *
 * @author Kangqi Fu
 */
public class VariableCopier {
    private VariableCopier(){
    }
    
    public static void copyToZerothRKStep(Cell[] cells){
        for (int i = Info.NUM_GHOST_CELLS; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS; i++){
            cells[i].u[0] = cells[i].u[Info.NUM_RK_STEPS];
        }
//        for (Cell cell : cells) {
//            cell.u[0] = cell.u[Info.NUM_RK_STEPS];
//        }
    }
}
