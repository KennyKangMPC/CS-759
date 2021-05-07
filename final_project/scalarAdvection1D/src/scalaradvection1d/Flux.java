package scalaradvection1d;

/**
 *
 * @author kenny
 */
public class Flux {
    private Flux(){
    }
    
    public static void calculateFlux(Cell[] cells, int rkStep){
        for (Cell cell : cells){
            cell.totalFlux = 0.0;
        }
        
        for (int interfaceIndex = Info.NUM_GHOST_CELLS; interfaceIndex < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS + 1; interfaceIndex++){
            double leftValue = cells[interfaceIndex - 1].uEast;
            double rightValue = cells[interfaceIndex].uWest;
            
            // Lax-Friedrichs scheme
            double flux = 0.5 * Info.ADVECTION_VEL * (leftValue + rightValue)
                    - 0.5 * Math.abs(Info.ADVECTION_VEL) * (rightValue - leftValue);
            cells[interfaceIndex - 1].totalFlux -= flux;
            cells[interfaceIndex].totalFlux += flux;
        }
    }
}
