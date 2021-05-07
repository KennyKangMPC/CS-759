package scalaradvection1d;

/**
 *
 * @author kenny
 */
public class GhostCellUpdater {
    private GhostCellUpdater(){
    }
    
    public static void updateGhostCells(Cell[] cells, int rkStep){
        for (int ghostCell = 0; ghostCell < Info.NUM_GHOST_CELLS; ghostCell++){
            cells[ghostCell].u[rkStep] = cells[Info.NUM_X_CELLS + ghostCell].u[rkStep];
            cells[Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS + ghostCell].u[rkStep]
                    = cells[Info.NUM_GHOST_CELLS + ghostCell].u[rkStep];
            
        }
    }
}
