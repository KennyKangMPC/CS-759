package shallowwater;

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
                for (int var = 0; var < Cell.NUM_VARS; var++)
                cell.totalFlux[var] = 0.0;
            }
        }
        
        //Assuming dx. dy are constant for the complete mesh
        double dx = cells[0][0].dx;
        double dy = cells[0][0].dy;
        
        // Flux calculatino for all vertical edges
        for (int j = Info.NUM_GHOST_CELLS; j < Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS; j ++) {
            for (int verInterfaceIdx = Info.NUM_GHOST_CELLS; verInterfaceIdx < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS + 1; verInterfaceIdx++){
                double[] leftValue = cells[verInterfaceIdx - 1][j].uEast;
                double[] rightValue = cells[verInterfaceIdx][j].uWest;
                
                double gravity = 9.8;
                double hL = leftValue[0];
                double uL = leftValue[1]/hL;
                double vL = leftValue[2]/hL;
                double cL = Math.sqrt(gravity * hL);
                double maxSpeedL = Math.abs(uL) + cL;
                
                double hR = rightValue[0];
                double uR = rightValue[1]/hR;
                double vR = rightValue[2]/hR;
                double cR = Math.sqrt(gravity * hR);
                double maxSpeedR = Math.abs(uR) + cR;
                
                double maxSpeed = Math.max(maxSpeedL, maxSpeedR);
                
                
                // Lax-Friedrichs scheme
                double[] FL = F(leftValue);
                double[] FR = F(rightValue);
                for (int var = 0; var < Cell.NUM_VARS; var++){
                    double flux = 0.5 * (FL[var] + FR[var]) - 0.5 * maxSpeed * (rightValue[var] - leftValue[var]);
                    flux *= dy;
                    cells[verInterfaceIdx-1][j].totalFlux[var] -= flux;
                    cells[verInterfaceIdx][j].totalFlux[var] += flux;
                }
            }
        }
        
        // Flux calculation for all horizontal edges
        for (int i = Info.NUM_GHOST_CELLS; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS; i ++) {
            for (int horInterfaceIdx = Info.NUM_GHOST_CELLS; horInterfaceIdx < Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS + 1; horInterfaceIdx++){
                double[] leftValue = cells[i][horInterfaceIdx - 1].uNorth;
                double[] rightValue = cells[i][horInterfaceIdx].uSouth;
                
                double gravity = 9.8;
                double hL = leftValue[0];
                //double uL = leftValue[1]/hL;
                double vL = leftValue[2]/hL;
                double cL = Math.sqrt(gravity * hL);
                double maxSpeedL = Math.abs(vL) + cL;
                
                double hR = rightValue[0];
                //double uR = rightValue[1]/hR;
                double vR = rightValue[2]/hR;
                double cR = Math.sqrt(gravity * hR);
                double maxSpeedR = Math.abs(vR) + cR;
                
                double maxSpeed = Math.max(maxSpeedL, maxSpeedR);
                
                double[] GL = G(leftValue);
                double[] GR = G(rightValue);
                for (int var = 0; var < Cell.NUM_VARS; var++){
                    
                    double flux = 0.5 * (GL[var]+GR[var]) - 0.5 * maxSpeed * (rightValue[var] - leftValue[var]);
                    flux *= dx;
                    cells[i][horInterfaceIdx - 1].totalFlux[var] -= flux;
                    cells[i][horInterfaceIdx].totalFlux[var] += flux;
                }
            }
        }

    }
    
    private static double[] F(double[] U) {
        double gravity = 9.81;
        double h = U[0];
        double u = U[1] / h;
        double v = U[2] / h;
        
        return new double[]{h*u, h*u*u+0.5*gravity*h*h, h*u*v};
    }
    
    private static double[] G(double[] U) {
        double gravity = 9.81;
        double h = U[0];
        double u = U[1] / h;
        double v = U[2] / h;
        
        return new double[]{h*v, h*u*v, h*v*v+0.5*gravity*h*h};
    }
}
