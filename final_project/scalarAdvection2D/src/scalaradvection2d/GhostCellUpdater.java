package scalaradvection2d;

/**
 *
 * @author kenny
 */
public class GhostCellUpdater {
    private GhostCellUpdater(){
    }
    
    public static void updateGhostCells(Cell[][] cells, int rkStep){
        // left and right side of the domain
        for (int j = 0; j < cells[0].length; j++){
            for (int ghostCell = 0; ghostCell < Info.NUM_GHOST_CELLS; ghostCell++){
                cells[ghostCell][j].u[rkStep] = cells[Info.NUM_X_CELLS + ghostCell][j].u[rkStep];
                cells[Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS + ghostCell][j].u[rkStep]
                    = cells[Info.NUM_GHOST_CELLS + ghostCell][j].u[rkStep];
            }
        }
        
        // bottom and top cells
        for (int i = 0; i < cells.length; i++){
            for (int ghostCell = 0; ghostCell < Info.NUM_GHOST_CELLS; ghostCell++){
                cells[i][ghostCell].u[rkStep] = cells[i][Info.NUM_Y_CELLS + ghostCell].u[rkStep];
                cells[i][Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS + ghostCell].u[rkStep]
                    = cells[i][Info.NUM_GHOST_CELLS + ghostCell].u[rkStep];
            }
        }
        

    }
}
