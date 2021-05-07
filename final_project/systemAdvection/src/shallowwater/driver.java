package shallowwater;

/**
 *
 * @author Kangqi Fu
 */
public class driver {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        
        System.out.println(Info.infoToString()); // print out the field information
        
        Cell[][] cells = new Cell[Info.NUM_X_CELLS + 2 * Info.NUM_GHOST_CELLS][Info.NUM_Y_CELLS + 2 * Info.NUM_GHOST_CELLS];
        double dx = (Info.MAX_X - Info.MIN_X) / Info.NUM_X_CELLS;
        double dy = (Info.MAX_Y - Info.MIN_Y) / Info.NUM_Y_CELLS;
        
        for (int i = 0; i < cells.length; i++){
            for (int j = 0; j < cells[0].length; j++){
                cells[i][j] = new Cell();
                cells[i][j].dx = dx;
                cells[i][j].dy = dy;
                cells[i][j].cx = Info.MIN_X + dx * (i + 0.5 - Info.NUM_GHOST_CELLS);
                cells[i][j].cy = Info.MIN_Y + dy * (j + 0.5 - Info.NUM_GHOST_CELLS);
            }
        }
        
        SchemeInitializer.initializeSolution(cells);
        double time = 0.0;
        SolutionFileWriter.writeSolutionFile(cells, time, Info.CFL);
        boolean lastTimeStep = false;
        
        long startTime = System.currentTimeMillis();
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
            System.out.printf("%d : Time = %f\n", timeIter, time);
            if (lastTimeStep){
                break;
            }
            
            if ((timeIter + 1) % 10 == 0)
                 SolutionFileWriter.writeSolutionFile(cells, time, Info.CFL);
        }
        long endTime = System.currentTimeMillis();
        System.out.println("Compute Time = " + (endTime - startTime));
        //System.out.println(new ErrorCalculator(cells)); // need to update this later
        
        System.out.println(Info.infoToString()); // print out the field information
        SolutionFileWriter.writeSolutionFile(cells, time, Info.CFL);
    }
    
}
