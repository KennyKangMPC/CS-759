package shallowwater;

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
                // recall this is periodic boundary condition
                for (int var = 0; var < Cell.NUM_VARS; var++){
                    cells[ghostCell][j].U[rkStep][var]
                        = Math.pow(-1, var) * cells[2 * Info.NUM_GHOST_CELLS - ghostCell - 1][j].U[rkStep][var];
                    cells[Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS + ghostCell][j].U[rkStep][var]
                        = Math.pow(-1, var) * cells[Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS - ghostCell - 1][j].U[rkStep][var];
                }
            }
        }
        
        // bottom and top cells
        for (Cell[] cell : cells) {
            for (int ghostCell = 0; ghostCell < Info.NUM_GHOST_CELLS; ghostCell++) {
                int var = 0;
                cell[ghostCell].U[rkStep][var] = cell[2 * Info.NUM_GHOST_CELLS - ghostCell - 1].U[rkStep][var];
                cell[Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS + ghostCell].U[rkStep][var] = cell[Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS - ghostCell - 1].U[rkStep][var];
                var = 1;
                cell[ghostCell].U[rkStep][var] = cell[2 * Info.NUM_GHOST_CELLS - ghostCell - 1].U[rkStep][var];
                cell[Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS + ghostCell].U[rkStep][var] = cell[Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS - ghostCell - 1].U[rkStep][var];
                var = 2;
                cell[ghostCell].U[rkStep][var] = -cell[2 * Info.NUM_GHOST_CELLS - ghostCell - 1].U[rkStep][var];
                cell[Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS + ghostCell].U[rkStep][var] = -cell[Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS - ghostCell - 1].U[rkStep][var];
            }
        }
//        // left and right side of the domain
//        for (int j = 0; j < cells[0].length; j++){
//            for (int ghostCell = 0; ghostCell < Info.NUM_GHOST_CELLS; ghostCell++){
//                
//                for (int var = 0; var < Cell.NUM_VARS; var++){
//                    cells[ghostCell][j].U[rkStep][var] = cells[Info.NUM_X_CELLS + ghostCell][j].U[rkStep][var];
//                    cells[Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS + ghostCell][j].U[rkStep][var]
//                        = cells[Info.NUM_GHOST_CELLS + ghostCell][j].U[rkStep][var];
//                }
//                
//            }
//        }
//        
//        // bottom and top cells
//        for (Cell[] cell : cells) {
//            for (int ghostCell = 0; ghostCell < Info.NUM_GHOST_CELLS; ghostCell++) {
//                for (int var = 0; var < Cell.NUM_VARS; var++) {
//                    cell[ghostCell].U[rkStep][var] = cell[Info.NUM_Y_CELLS + ghostCell].U[rkStep][var];
//                    cell[Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS + ghostCell].U[rkStep][var] = cell[Info.NUM_GHOST_CELLS + ghostCell].U[rkStep][var];
//                }
//            }
//        }
        
        
    }
}
