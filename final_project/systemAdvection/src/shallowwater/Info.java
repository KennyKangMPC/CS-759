package shallowwater;

import java.lang.reflect.Field;
import java.util.Arrays;
import java.util.stream.Collectors;

/**
 *
 * @author kenny
 */
public class Info {
    public static final int NUM_X_CELLS = 80;
    public static final int NUM_Y_CELLS = 80;
    
    public static final int NUM_GHOST_CELLS = 4;
    public static final int NUM_RK_STEPS = 3; // can be 1, 2, or 3 here
    public static final int MAX_TIME_ITER = 100000;
    
    public static final double STOPPING_TIME = 5.0;
    public static final double MIN_X = -1.0;
    public static final double MAX_X = 1.0;
    public static final double MIN_Y = -1.0;
    public static final double MAX_Y = 1.0;
    
    
    public static final double ADVECTION_VEL_X = 1.0;
    public static final double ADVECTION_VEL_Y = 1.0;
    public static final double CFL = 0.5;
    public static final ReconstructionTypes RECONST_TYPE = ReconstructionTypes.FIRST_ORDER;
    
    public static final String BASEPATH = "./output/";
    
    // For printing out the fields of the info file (uses Java )
    
    static String infoToString(){
        Field[] fields = Info.class.getFields();
        String str = Arrays.stream(fields)
                .map(field -> formatField(field))
                .collect(Collectors.joining("\n"));
        return str;
    }
    
    private static String formatField(Field field){
        String str = "";
        try {
            str = field.getName() + " : " + field.get(null);
        } catch (IllegalAccessException ex){
        }
        return str;
    }
    
}
