package shallowwater;

/**
 *
 * @author Kangqi Fu
 */
public class TimeIntegration {
    
    private TimeIntegration(){
    }
    
    public static void updateCellAverages(Cell[][] cells, int rkStep, double dt){
        // assuming dx, dy are constants
        double area = cells[0][0].dx * cells[0][0].dy;
        switch (Info.NUM_RK_STEPS){
            case 1: // 1 Step RK method
                for (int i = Info.NUM_GHOST_CELLS; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS; i++){
                    for (int j = Info.NUM_GHOST_CELLS; j < Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS; j++)
                        for (int var = 0; var < Cell.NUM_VARS; var++){
                            cells[i][j].U[rkStep + 1][var] = cells[i][j].U[rkStep][var] + dt / area * cells[i][j].totalFlux[var];
                        }
                }
                break;
            case 2: // 2 Step RK method
                for (int i = Info.NUM_GHOST_CELLS; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS; i++) {
                    for (int j = Info.NUM_GHOST_CELLS; j < Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS; j++){
                        for (int var = 0; var < Cell.NUM_VARS; var++){
                            switch (rkStep){
                                case 0: //step 0
                                        cells[i][j].U[rkStep + 1][var] 
                                            = cells[i][j].U[rkStep][var] 
                                            + dt / area * cells[i][j].totalFlux[var];
                                    break;
                                case 1: // step 1
                                        cells[i][j].U[rkStep + 1][var] 
                                            = 0.5 * (cells[i][j].U[rkStep - 1][var] 
                                            + cells[i][j].U[rkStep][var] 
                                            + dt / area * cells[i][j].totalFlux[var]);
                                    break;
                            default:
                                throw new IllegalStateException("Shouldn't reach here! Errors exist!");
                            }
                        }
                    }
                }
                break;
                
            case 3: // 3 step RK method
                for (int i = Info.NUM_GHOST_CELLS; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS; i++) {
                    for (int j = Info.NUM_GHOST_CELLS; j < Info.NUM_Y_CELLS + Info.NUM_GHOST_CELLS; j++){
                        for (int var = 0; var < Cell.NUM_VARS; var++){
                            switch (rkStep) {
                                case 0:
                                    cells[i][j].U[rkStep + 1][var]
                                            = cells[i][j].U[rkStep][var]
                                            + dt / area * cells[i][j].totalFlux[var];
                                    break;
                                case 1:
                                    cells[i][j].U[rkStep + 1][var]
                                            = 3.0/4.0 * cells[i][j].U[rkStep - 1][var]
                                            + 1.0/4.0 * cells[i][j].U[rkStep][var]
                                            + 1.0/4.0 * dt/area * cells[i][j].totalFlux[var];
                                    break;
                                case 2:
                                    cells[i][j].U[rkStep + 1][var]
                                            = 1.0/3.0 * cells[i][j].U[rkStep - 2][var]
                                            + 2.0/3.0 * cells[i][j].U[rkStep][var]
                                            + 2.0/3.0 * dt/area * cells[i][j].totalFlux[var];
                                    break;
                                default:
                                    throw new IllegalStateException("Shouldn't reach here! Errors exist!");
                            }
                        }
                    }
                }
                break; //TODO: Add 4th order
            default:
                System.out.println("The " + Info.NUM_RK_STEPS + "step RK method has not been defined yet!");
                throw new RuntimeException("The " + Info.NUM_RK_STEPS + "step RK method has not been defined yet!");
        }
        
    }
}
