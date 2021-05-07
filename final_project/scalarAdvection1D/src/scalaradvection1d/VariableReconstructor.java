package scalaradvection1d;

import weno.ReconstructedValues;
import weno.WENO;

/**
 *
 * @author kenny
 */
public class VariableReconstructor {
    private VariableReconstructor(){
    }
    
    public static void reconstructVariables(Cell[] cells, int rkStep){
        switch(Info.RECONST_TYPE) {
            case FIRST_ORDER:
                for (int i = Info.NUM_GHOST_CELLS-1 ; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS + 1; i++) {
                    cells[i].uWest = cells[i].u[rkStep];
                    cells[i].uEast = cells[i].u[rkStep];
                }
                break;
            case BEAM_WARMING: // numerical oscillations, need limiters, diffustions are reduced 
                for (int i = Info.NUM_GHOST_CELLS - 1; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS + 1; i++) {
                    double du_dx = (cells[i].u[rkStep] - cells[i-1].u[rkStep]) / cells[i].dx;
                    cells[i].uWest = cells[i].u[rkStep] - du_dx * cells[i].dx/2.0;
                    cells[i].uEast = cells[i].u[rkStep] + du_dx * cells[i].dx/2.0;
                }
                break;
            case LAX_WENDROFF: // more osciallation
                for (int i = Info.NUM_GHOST_CELLS - 1; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS + 1; i++) {
                    double du_dx = (cells[i + 1].u[rkStep] - cells[i].u[rkStep]) / cells[i].dx;
                    cells[i].uWest = cells[i].u[rkStep] - du_dx * cells[i].dx/2.0;
                    cells[i].uEast = cells[i].u[rkStep] + du_dx * cells[i].dx/2.0;
                }
                break;
            case FROMM: // less osciallation; different CFL value
                for (int i = Info.NUM_GHOST_CELLS - 1; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS + 1; i++) {
                    double du_dx = (cells[i + 1].u[rkStep] - cells[i - 1].u[rkStep]) / 2.0 / cells[i].dx;
                    cells[i].uWest = cells[i].u[rkStep] - du_dx * cells[i].dx/2.0;
                    cells[i].uEast = cells[i].u[rkStep] + du_dx * cells[i].dx/2.0;
                }
                break;
            case LIMITED_LW:
                for (int i = Info.NUM_GHOST_CELLS - 1; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS + 1; i++) {
                    
                    double r = (cells[i].u[rkStep] - cells[i-1].u[rkStep] + 1e-6)
                            / (cells[i+1].u[rkStep]-cells[i].u[rkStep] + 1e-6);
                    r = Math.max(0, r);
                    // use Albada limiter. there are many other possible limiters can be used: 
                    double phi = (r * r + r) / (r * r + 1.0);
                    double du_dx = (cells[i + 1].u[rkStep] - cells[i].u[rkStep]) / cells[i].dx;
                    cells[i].uWest = cells[i].u[rkStep] - phi * du_dx * cells[i].dx/2.0;
                    cells[i].uEast = cells[i].u[rkStep] + phi * du_dx * cells[i].dx/2.0;
                }
                break;
                
            case WENO:
                for (int i = Info.NUM_GHOST_CELLS - 1; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS + 1; i++){
                    int k = 3; // k = 2 gives WENO of order 3, k = 3, gives order of 5
                    int numNeigh = k - 1;
                    double[] averageValues = new double[2 * k - 1];
                    for (int j = -numNeigh; j <= numNeigh; j++){
                        averageValues[j + numNeigh] = cells[i + j].u[rkStep];
                    }
                    ReconstructedValues values = new WENO().getWenoReconstructedValues(averageValues, k);
                    cells[i].uWest = values.value_imh; // i minus h
                    cells[i].uEast = values.value_iph; // i minus h
                }
                break;
            default:
                throw new UnsupportedOperationException("This reconstruction type is not defined yet");
        }
        
    }
}
