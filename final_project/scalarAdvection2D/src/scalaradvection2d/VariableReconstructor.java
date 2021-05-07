package scalaradvection2d;

import java.util.stream.IntStream;
import weno.ReconstructedValues;
import weno.WENO;

/**
 *
 * @author kenny
 */
public class VariableReconstructor {
    private VariableReconstructor(){
    }
    
    public static void reconstructVariables(Cell[][] cells, int rkStep){
        switch(Info.RECONST_TYPE) {
            case FIRST_ORDER:
                for (int i = Info.NUM_GHOST_CELLS-1 ; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS + 1; i++) {
                    for (int j = Info.NUM_GHOST_CELLS-1 ; j < Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS + 1; j++) {
                        cells[i][j].uWest = cells[i][j].u[rkStep];
                        cells[i][j].uEast = cells[i][j].u[rkStep];
                        cells[i][j].uSouth = cells[i][j].u[rkStep];
                        cells[i][j].uNorth = cells[i][j].u[rkStep];
                    }
                }
                break;
            case BEAM_WARMING: // numerical oscillations, need limiters, diffustions are reduced 
                
                for (int i = Info.NUM_GHOST_CELLS - 1; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS + 1; i++) {
                    for (int j = Info.NUM_GHOST_CELLS - 1; j < Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS + 1; j++) {
                        double du_dx = (cells[i][j].u[rkStep] - cells[i-1][j].u[rkStep]) / cells[i][j].dx;
                        double du_dy = (cells[i][j].u[rkStep] - cells[i][j-1].u[rkStep]) / cells[i][j].dy;
                        cells[i][j].uWest = cells[i][j].u[rkStep] - du_dx * cells[i][j].dx/2.0;
                        cells[i][j].uEast = cells[i][j].u[rkStep] + du_dx * cells[i][j].dx/2.0;
                        
                        cells[i][j].uSouth = cells[i][j].u[rkStep] - du_dy * cells[i][j].dy/2.0;
                        cells[i][j].uNorth = cells[i][j].u[rkStep] + du_dy * cells[i][j].dy/2.0;
                    }
                }
                break;
            case LAX_WENDROFF: // more osciallation
                for (int i = Info.NUM_GHOST_CELLS - 1; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS + 1; i++) {
                    for (int j = Info.NUM_GHOST_CELLS - 1; j < Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS + 1; j++) {
                        double du_dx = (cells[i+1][j].u[rkStep] - cells[i][j].u[rkStep]) / cells[i][j].dx;
                        double du_dy = (cells[i][j+1].u[rkStep] - cells[i][j].u[rkStep]) / cells[i][j].dy;
                        cells[i][j].uWest = cells[i][j].u[rkStep] - du_dx * cells[i][j].dx/2.0;
                        cells[i][j].uEast = cells[i][j].u[rkStep] + du_dx * cells[i][j].dx/2.0;
                        
                        cells[i][j].uSouth = cells[i][j].u[rkStep] - du_dy * cells[i][j].dy/2.0;
                        cells[i][j].uNorth = cells[i][j].u[rkStep] + du_dy * cells[i][j].dy/2.0;
                    }
                }
                break;
            case FROMM: // less osciallation; different CFL value
                for (int i = Info.NUM_GHOST_CELLS - 1; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS + 1; i++) {
                    for (int j = Info.NUM_GHOST_CELLS - 1; j < Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS + 1; j++) {
                        double du_dx = (cells[i+1][j].u[rkStep] - cells[i-1][j].u[rkStep]) / 2.0 / cells[i][j].dx;
                        double du_dy = (cells[i][j+1].u[rkStep] - cells[i][j-1].u[rkStep]) / 2.0 / cells[i][j].dy;
                        cells[i][j].uWest = cells[i][j].u[rkStep] - du_dx * cells[i][j].dx/2.0;
                        cells[i][j].uEast = cells[i][j].u[rkStep] + du_dx * cells[i][j].dx/2.0;
                        
                        cells[i][j].uSouth = cells[i][j].u[rkStep] - du_dy * cells[i][j].dy/2.0;
                        cells[i][j].uNorth = cells[i][j].u[rkStep] + du_dy * cells[i][j].dy/2.0;
                    }
                }
                break;
            case LIMITED_LW:
                for (int i = Info.NUM_GHOST_CELLS - 1; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS + 1; i++) {
                    for (int j = Info.NUM_GHOST_CELLS - 1; j < Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS + 1; j++) {
                        
                        double numerator = cells[i][j].u[rkStep] - cells[i-1][j].u[rkStep];
                        double denominator = cells[i+1][j].u[rkStep] - cells[i][j].u[rkStep];
                        
                        double r_x = numerator / denominator;
                        
                        if (Math.abs(denominator) < 1e-15){
                            r_x = 0.0;
                        }
                        r_x = Math.max(0, r_x);
                        double phi_x = (r_x * r_x + r_x) / (r_x * r_x + 1.0);
                        
                        numerator = cells[i][j].u[rkStep] - cells[i][j-1].u[rkStep];
                        denominator = cells[i][j+1].u[rkStep] - cells[i][j].u[rkStep];
                        
                        double r_y = numerator / denominator;
                        if (Math.abs(denominator) < 1e-15) {
                            r_y = 0.0;
                        }
                        r_y = Math.max(0, r_y);
                        double phi_y = (r_y * r_y + r_y) / (r_y * r_y + 1.0);
                        
                        double du_dx = (cells[i+1][j].u[rkStep] - cells[i][j].u[rkStep]) / cells[i][j].dx;
                        double du_dy = (cells[i][j+1].u[rkStep] - cells[i][j].u[rkStep]) / cells[i][j].dy;
                        
                        cells[i][j].uWest = cells[i][j].u[rkStep] - phi_x * du_dx * cells[i][j].dx/2.0;
                        cells[i][j].uEast = cells[i][j].u[rkStep] + phi_x * du_dx * cells[i][j].dx/2.0;
                        
                        cells[i][j].uSouth = cells[i][j].u[rkStep] - phi_y * du_dy * cells[i][j].dy/2.0;
                        cells[i][j].uNorth = cells[i][j].u[rkStep] + phi_y * du_dy * cells[i][j].dy/2.0;
                    }
                }
                
                break;
            case LIMITED_BW:
                for (int i = Info.NUM_GHOST_CELLS - 1; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS + 1; i++) {
                    for (int j = Info.NUM_GHOST_CELLS - 1; j < Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS + 1; j++) {
                        
                        double numerator = cells[i+1][j].u[rkStep] - cells[i][j].u[rkStep];
                        double denominator = cells[i][j].u[rkStep] - cells[i-1][j].u[rkStep];
                        
                        double r_x = numerator / denominator;
                        
                        if (Math.abs(denominator) < 1e-15){
                            r_x = 0.0;
                        }
                        r_x = Math.max(0, r_x);
                        double phi_x = (r_x * r_x + r_x) / (r_x * r_x + 1.0);
                        
                        numerator = cells[i][j+1].u[rkStep] - cells[i][j].u[rkStep];
                        denominator = cells[i][j].u[rkStep] - cells[i][j-1].u[rkStep];
                        
                        double r_y = numerator / denominator;
                        if (Math.abs(denominator) < 1e-15) {
                            r_y = 0.0;
                        }
                        r_y = Math.max(0, r_y);
                        double phi_y = (r_y * r_y + r_y) / (r_y * r_y + 1.0);
                        
                        double du_dx = (cells[i][j].u[rkStep] - cells[i-1][j].u[rkStep]) / cells[i][j].dx;
                        double du_dy = (cells[i][j].u[rkStep] - cells[i][j-1].u[rkStep]) / cells[i][j].dy;
                        
                        cells[i][j].uWest = cells[i][j].u[rkStep] - phi_x * du_dx * cells[i][j].dx/2.0;
                        cells[i][j].uEast = cells[i][j].u[rkStep] + phi_x * du_dx * cells[i][j].dx/2.0;
                        
                        cells[i][j].uSouth = cells[i][j].u[rkStep] - phi_y * du_dy * cells[i][j].dy/2.0;
                        cells[i][j].uNorth = cells[i][j].u[rkStep] + phi_y * du_dy * cells[i][j].dy/2.0;
                    }
                }
                break;
            case WENO: // will generalized to // This way makes the code much faster !!! I have 4 processors //weno paper - Weighted essentially non oscillatory
                IntStream.range(Info.NUM_GHOST_CELLS - 1, Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS + 1)
                        .parallel()
                        .forEach(i -> {
                            for (int j = Info.NUM_GHOST_CELLS - 1; j < Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS + 1; j++) {
                                int k = 3; // k = 2 gives WENO of order 3, k = 3, gives order of 5, k = 4 gives order of 7
                                int numNeigh = k - 1;
                                double[] averageValues = new double[2 * k - 1];
                                for (int n = -numNeigh; n <= numNeigh; n++){
                                    averageValues[n + numNeigh] = cells[i + n][j].u[rkStep];
                                }
                                ReconstructedValues values = new WENO().getWenoReconstructedValues(averageValues, k);
                                cells[i][j].uWest = values.value_imh; // i minus h
                                cells[i][j].uEast = values.value_iph; // i minus h
                        
                                for (int n = -numNeigh; n <= numNeigh; n++){
                                    averageValues[n + numNeigh] = cells[i][j+n].u[rkStep];
                                }
                                values = new WENO().getWenoReconstructedValues(averageValues, k);
                                cells[i][j].uSouth = values.value_imh; // i minus h
                                cells[i][j].uNorth = values.value_iph; // i plus h
                            }
                        });
                break;
                        
                  // below are the       
//                for (int i = Info.NUM_GHOST_CELLS - 1; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS + 1; i++){
//                    for (int j = Info.NUM_GHOST_CELLS - 1; j < Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS + 1; j++) {
//                        int k = 3; // k = 2 gives WENO of order 3, k = 3, gives order of 5, k = 4 gives order of 7
//                        int numNeigh = k - 1;
//                        double[] averageValues = new double[2 * k - 1];
//                        for (int n = -numNeigh; n <= numNeigh; n++){
//                            averageValues[n + numNeigh] = cells[i + n][j].u[rkStep];
//                        }
//                        ReconstructedValues values = new WENO().getWenoReconstructedValues(averageValues, k);
//                        cells[i][j].uWest = values.value_imh; // i minus h
//                        cells[i][j].uEast = values.value_iph; // i minus h
//                        
//                        for (int n = -numNeigh; n <= numNeigh; n++){
//                            averageValues[n + numNeigh] = cells[i][j+n].u[rkStep];
//                        }
//                        values = new WENO().getWenoReconstructedValues(averageValues, k);
//                        cells[i][j].uSouth = values.value_imh; // i minus h
//                        cells[i][j].uNorth = values.value_iph; // i plus h
//                    }
//                }
//                break;
            case SDWLS: // solution dependent weighted least square 
//                J. C. Mandal and S. P. Rao, High resolution nite volume computations on unstructured
//                grids using solution dependent weighted least squares gradients, Computers and Fluids, 
//                vol. 44, pp. 2331, May 2011
                
                // one type of WENO?
                final double EPS = 1e-15;
                for (int i = Info.NUM_GHOST_CELLS - 1; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS + 1; i++) {
                    for (int j = Info.NUM_GHOST_CELLS - 1; j < Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS + 1; j++) {
                        double u_ij = cells[i][j].u[rkStep];
                        double u_ip1j = cells[i + 1][j].u[rkStep]; //u_i_plus_1
                        double u_im1j = cells[i - 1][j].u[rkStep]; // u_i_minus_1
                        double u_ijp1 = cells[i][j + 1].u[rkStep]; 
                        double u_ijm1 = cells[i][j - 1].u[rkStep];

                        double dx = cells[i][j].dx;
                        double dy = cells[i][j].dy;
                        
                        double du1 = u_ip1j - u_ij;
                        double du2 = u_im1j - u_ij;
                        double w1 = 1.0 / (du1 * du1 + EPS);
                        double w2 = 1.0 / (du2 * du2 + EPS);
                        
                        double du_dx = (w1 * du1 - w2 * du2) / (w1 + w2) / dx;
                        cells[i][j].uWest = cells[i][j].u[rkStep] - du_dx * dx / 2.0;
                        cells[i][j].uEast = cells[i][j].u[rkStep] + du_dx * dx / 2.0;
                        
                        du1 = u_ijp1 - u_ij;
                        du2 = u_ijm1 - u_ij;
                        w1 = 1.0 / (du1 * du1 + EPS);
                        w2 = 1.0 / (du2 * du2 + EPS);
                        double du_dy = (w1 * du1 - w2 * du2) / (w1 + w2) / dy;
                        cells[i][j].uSouth = cells[i][j].u[rkStep] - du_dy * dy / 2.0;
                        cells[i][j].uNorth = cells[i][j].u[rkStep] + du_dy * dy / 2.0;
                    }
                }
                break;
            default:
                throw new UnsupportedOperationException("This reconstruction type is not defined yet");
        }
        
    }
}
