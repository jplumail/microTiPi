/**
 *
 */
package microTiPi.microUtils;

import mitiv.array.ArrayUtils;
import mitiv.base.Shape;
import mitiv.utils.Imager;
import mitiv.utils.TiPiHook;

/**
 * @author ferreol
 *
 */
public class DEMICSHook implements TiPiHook{

    private Imager curImager;
    private Shape outShape;
    private boolean debug;

    /**
     *
     */
    public DEMICSHook(Imager imager, Shape outShape, boolean debug) {
        this.curImager = imager;
        this.outShape = outShape;
        this.debug = debug;
    }
    /* (non-Javadoc)
     * @see mitiv.utils.TiPiHook#run(java.lang.Object, int)
     */
    @Override
    public void run(Object caller, int iter) {
        // TODO Auto-generated method stub
        curImager.show(ArrayUtils.crop(((MicroDeconvolution) caller).solver.getSolution(),outShape) ,"Current mu="+((MicroDeconvolution) caller).solver.getRegularizationLevel() +"it:"+((MicroDeconvolution) caller).solver.getIterations());

        if (debug){
            System.out.println("Cost "+((MicroDeconvolution) caller).solver.getCost() );
        }
    }


}
