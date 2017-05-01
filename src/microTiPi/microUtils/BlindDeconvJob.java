/**
 *
 */
package microTiPi.microUtils;

import microTiPi.microscopy.PSF_Estimation;
import mitiv.array.ArrayUtils;
import mitiv.array.ShapedArray;

/**
 * @author ferreol
 *
 */
public class BlindDeconvJob {

    private int totalNbOfBlindDecLoop;
    private ShapedArray psfArray;
    private boolean debug=false;
    //   private ShapedArray objArray;
    private PSF_Estimation psfEstimation;
    private DeconvolutionJob deconvolver;
    private int[] parametersFlags;
    private boolean run =true;
    private int[] maxIter;

    /**
     *
     */
    public BlindDeconvJob(int totalNbOfBlindDecLoop,int[] parametersFlags,int[] maxIter, PSF_Estimation psfEstimation ,DeconvolutionJob deconvolver, boolean debug ) {
        this.totalNbOfBlindDecLoop = totalNbOfBlindDecLoop;
        this.parametersFlags = parametersFlags;
        this.maxIter = maxIter;
        this.psfEstimation = psfEstimation;
        this.deconvolver = deconvolver;
        this.debug = debug;

    }

    public void blindDeconv(ShapedArray objArray){
        run =true;
        for(int i = 0; i < totalNbOfBlindDecLoop; i++) {
            psfArray = ArrayUtils.roll(psfEstimation.getPupil().getPsf());
            psfEstimation.freeMem();

            deconvolver.updatePsf(psfArray);


            objArray = deconvolver.deconv(objArray);

            //Emergency stop
            if (!run) {
                return;
            }
            psfEstimation.setObj(objArray);

            for (int j = 0; j < parametersFlags.length; j++) {
                if (debug ) {
                    System.out.println("------------------");
                    System.out.println("  "+ j+" estimation");
                    System.out.println("------------------");
                }
                psfEstimation.setRelativeTolerance(0.);
                psfEstimation.setMaximumIterations(maxIter[j]);
                psfEstimation.fitPSF( parametersFlags[j]);
                //Emergency stop
                if (!run) {
                    return;
                }
            }
        }

    }

    /**
     * Emergency stop
     */
    public void abort(){
        System.out.println("abort");
        run  = false;
        deconvolver.abort();
        psfEstimation.abort();
    }

}
