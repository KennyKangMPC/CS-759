package scalaradvection1d;

/**
 *
 * @author Kangqi Fu
 */
public class TimeIntegration {
    
    private TimeIntegration(){
    }
    
    public static void updateCellAverages(Cell[] cells, int rkStep, double dt){
        switch (Info.NUM_RK_STEPS){
            case 1: // 1 Step RK method
                for (int i = Info.NUM_GHOST_CELLS; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS; i++){
                    cells[i].u[rkStep + 1] = cells[i].u[rkStep] + dt / cells[i].dx * cells[i].totalFlux;
                }
                break;
            case 2: // 2 Step RK method
                for (int i = Info.NUM_GHOST_CELLS; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS; i++) {
                    switch (rkStep){
                        case 0: //step 0
                            cells[i].u[rkStep + 1] 
                                    = cells[i].u[rkStep] 
                                    + dt / cells[i].dx * cells[i].totalFlux;
                            break;
                        case 1: // step 1
                            cells[i].u[rkStep + 1] 
                                    = 0.5 * (cells[i].u[rkStep - 1] 
                                    + cells[i].u[rkStep] 
                                    + dt / cells[i].dx * cells[i].totalFlux);
                            break;
                        default:
                            throw new IllegalStateException("Shouldn't reach here! Errors exist!");
                    }
                }
                break;
            case 3: // 3 step RK method
                for (int i = Info.NUM_GHOST_CELLS; i < Info.NUM_X_CELLS + Info.NUM_GHOST_CELLS; i++) {
                    switch (rkStep) {
                        case 0:
                            cells[i].u[rkStep + 1]
                                    = cells[i].u[rkStep]
                                    + dt / cells[i].dx * cells[i].totalFlux;
                            break;
                        case 1:
                            cells[i].u[rkStep + 1]
                                    = 3.0/4.0 * cells[i].u[rkStep - 1]
                                    + 1.0/4.0 * cells[i].u[rkStep]
                                    + 1.0/4.0 * dt/cells[i].dx * cells[i].totalFlux;
                            break;
                        case 2:
                            cells[i].u[rkStep + 1]
                                    = 1.0/3.0 * cells[i].u[rkStep - 2]
                                    + 2.0/3.0 * cells[i].u[rkStep]
                                    + 2.0/3.0 * dt/cells[i].dx * cells[i].totalFlux;
                            break;
                        default:
                            throw new IllegalStateException("Shouldn't reach here! Errors exist!");
                    }
                }
                break; //TODO: Add 4th order
            default:
                System.out.println("The " + Info.NUM_RK_STEPS + "step RK method has not been defined yet!");
                throw new RuntimeException("The " + Info.NUM_RK_STEPS + "step RK method has not been defined yet!");
        }
        
    }
}
