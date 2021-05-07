package scalaradvection2d;

/**
 *
 * @author Kangqi Fu
 */
public class SchemeInitializer {
    private SchemeInitializer(){
    }
    
    public static void initializeSolution(Cell[][] cells){
        int rkStep = 0;
        
        for (Cell[] cellRow : cells){
            for (Cell cell : cellRow){
                if (cell.cx > -0.25 && cell.cx < 0.25 && cell.cy > -0.25 && cell.cy < 0.25){
                    cell.u[rkStep] = 1.0;
                } else {
                    cell.u[rkStep] = 0.0;
                }
            }
            //  /Sine wave
////            double x = cell.cx;
////            cell.u[rkStep] = Math.sin(2.0 * Math.PI * x);
//            
//            // Square wave
//            if (cell.cx > -0.25 && cell.cx < 0.25) {
//                cell.u[rkStep] = 1.0;
//            } else {
//                cell.u[rkStep] = 0.0;
//            }
//            
//            // Gaussian
//            double x = cell.cx;
//            double sigma = Math.sqrt(0.05);
//            double mu = 0.0;
//            double a = 1.0 / (sigma * Math.sqrt(2.0 * Math.PI));
//            double b = mu;
//            double c = sigma;
//            cell.u[rkStep] = a * Math.exp(-((x-b) * (x-b))) / (2.0*c*c);
//            
//            // Combination of sub-funciton
//            double rx = (cell.cx - Info.MIN_X) / (Info.MAX_X - Info.MIN_X);
//            if (rx >= 0.1 && rx <= 0.2){
//                cell.u[rkStep] = Math.exp(-1111.111111 * Math.log(2.0) * Math.pow((rx - 0.15) * 2, 2.0));
//            } else if (rx >= 0.3 && rx <= 0.4){
//                cell.u[rkStep] = 1.0;
//            } else if (rx >= 0.5 && rx <= 0.6) {
//                cell.u[rkStep] = 1.0 - 10 * Math.abs((rx - 0.55) * 2);
//            } else if (rx >= 0.7 && rx <= 0.8) {
//                cell.u[rkStep] = 1.0 - 100 * Math.pow((rx - 0.75) * 2, 2.0);
//            } else {
//                cell.u[rkStep] = 0.0;
//            }
        }
        
    }
}
