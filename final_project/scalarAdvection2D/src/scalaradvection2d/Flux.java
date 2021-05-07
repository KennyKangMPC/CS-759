package scalaradvection2d;

/**
 *
 * @author kenny
 */
public class Flux {
    private Flux(){
    }
    
    public static void calculateFlux(Cell[][] cells, int rkStep){
        for (Cell[] cellRow : cells) {
            for (Cell cell : cellRow){
                cell.totalFlux = 0.0;
            }
        }
        
        double dx = cells[0][0].dx;
        double dy = cells[0][0].dy;
        
        // Flux calculatino for all vertical edges
        for (int j = Info.NUM_GHOST_CELLS; j < Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS; j ++) {
            for (int verInterfaceIdx = Info.NUM_GHOST_CELLS; verInterfaceIdx < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS + 1; verInterfaceIdx++){
                double leftValue = cells[verInterfaceIdx - 1][j].uEast;
                double rightValue = cells[verInterfaceIdx][j].uWest;
            
                // Lax-Friedrichs scheme
                
                double flux = 0.5 * Info.ADVECTION_VEL_X * (leftValue + rightValue)
                    - 0.5 * Math.abs(Info.ADVECTION_VEL_X) * (rightValue - leftValue);
                flux *= dy;
                cells[verInterfaceIdx - 1][j].totalFlux -= flux;
                cells[verInterfaceIdx][j].totalFlux += flux;
            }
        }
        
        // Flux calculation for all horizontal edges
        for (int i = Info.NUM_GHOST_CELLS; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS; i ++) {
            for (int horInterfaceIdx = Info.NUM_GHOST_CELLS; horInterfaceIdx < Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS + 1; horInterfaceIdx++){
                double leftValue = cells[i][horInterfaceIdx - 1].uNorth;
                double rightValue = cells[i][horInterfaceIdx].uSouth;
            
                // Lax-Friedrichs scheme
                double flux = 0.5 * Info.ADVECTION_VEL_Y * (leftValue + rightValue)
                    - 0.5 * Math.abs(Info.ADVECTION_VEL_Y) * (rightValue - leftValue);
                flux *= dx;
                cells[i][horInterfaceIdx - 1].totalFlux -= flux;
                cells[i][horInterfaceIdx].totalFlux += flux;
            }
        }

    }
}
