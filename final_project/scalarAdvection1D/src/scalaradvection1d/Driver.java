package scalaradvection1d;

/**
 * This is the main class of this library. 
 * @author Kangqi Fu
 */
public class Driver {
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        Cell[] cells = new Cell[Info.NUM_X_CELLS + 2 * Info.NUM_GHOST_CELLS];
        double dx = (Info.MAX_X - Info.MIN_X) / Info.NUM_X_CELLS;
        for (int i = 0; i < cells.length; i++){
            cells[i] = new Cell();
            cells[i].dx = dx;
            cells[i].cx = Info.MIN_X + dx * (i + 0.5 - Info.NUM_GHOST_CELLS);
        }
        
        SchemeInitializer.initializeSolution(cells);
        double time = 0.0;
        SolutionFileWriter.writeSolutionFile(cells, time, Info.CFL);
        boolean lastTimeStep = false;
        
        for (int timeIter = 0; timeIter < Info.MAX_TIME_ITER; timeIter++) {
            double dt = TimeStepCalculator.getTimeStep(cells);
            if (time + dt > Info.STOPPING_TIME) {
                dt = Info.STOPPING_TIME - time;
                lastTimeStep = true;
            }
            
            for (int rkStep = 0; rkStep < Info.NUM_RK_STEPS; rkStep++){
                GhostCellUpdater.updateGhostCells(cells, rkStep);
                VariableReconstructor.reconstructVariables(cells, rkStep);
                Flux.calculateFlux(cells, rkStep);
                TimeIntegration.updateCellAverages(cells, rkStep, dt);
            }
            
            VariableCopier.copyToZerothRKStep(cells);
            
            time += dt;
            if (lastTimeStep){
                break;
            }
            
            
            if ((timeIter + 1) % 5 == 0)
                 SolutionFileWriter.writeSolutionFile(cells, time, Info.CFL);
        }
        
        System.out.println("The numerical used here is " + Info.RECONST_TYPE);
        System.out.println(new ErrorCalculator(cells));
        SolutionFileWriter.writeSolutionFile(cells, time, Info.CFL);
    }
    
}
