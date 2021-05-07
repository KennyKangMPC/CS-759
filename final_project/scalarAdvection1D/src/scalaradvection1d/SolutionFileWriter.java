package scalaradvection1d;
import java.io.File;
import java.io.IOException;
import java.io.FileWriter;

/**
 *
 * @author Kangqi Fu
 */
public class SolutionFileWriter {
    static int fileCounter = 0;
    private SolutionFileWriter(){
    }
    
    static void writeSolutionFile(Cell[] cells, double time, double cfl){
        int rkStep = 0;
        File file = new File(Info.BASEPATH, String.format("Solution%04d.dat", fileCounter));
        try{
            try(FileWriter fileWriter = new FileWriter(file)){
                fileWriter.write("NumXCells:" + Info.NUM_X_CELLS + "\n");
                fileWriter.write("NumGhostCells:" + Info.NUM_GHOST_CELLS + "\n");
                fileWriter.write("time:" + time + "\n");
                fileWriter.write("CFL:" + cfl + "\n");
                for (Cell cell : cells){
                    fileWriter.write(String.format("%20.15f %20.15f\n", cell.cx, cell.u[rkStep]));
                }
            }
        } catch (IOException ex){
            System.out.println("Unable to write solution file");
            System.out.println(ex.getMessage());
        }
        
        fileCounter++;
    }
    
}
