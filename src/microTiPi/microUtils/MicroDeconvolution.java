/**
 *
 */
package microTiPi.microUtils;

import mitiv.array.ShapedArray;
import mitiv.base.Shape;
import mitiv.cost.EdgePreservingDeconvolution;
import mitiv.optim.OptimTask;

/**
 * @author ferreol
 *
 */
public class MicroDeconvolution {
    protected EdgePreservingDeconvolution solver;
    protected DEMICSHook iterHook,finalHook;
    protected boolean run=true;

    public MicroDeconvolution(ShapedArray dataArray, ShapedArray psfArray, ShapedArray wgtArray,
            Shape outputShape, double mu, double epsilon, double[] scale, boolean positivity,
            boolean single, int nbIterDeconv, DEMICSHook iterHook, DEMICSHook finalHook){

        solver = new EdgePreservingDeconvolution();

        solver.setForceSinglePrecision(single);

        solver.setRelativeTolerance(0.0);
        solver.setUseNewCode(false);
        solver.setObjectShape(outputShape);
        solver.setPSF(psfArray);
        solver.setData(dataArray);
        solver.setWeights(wgtArray);
        solver.setEdgeThreshold(epsilon);
        solver.setRegularizationLevel(mu);

        solver.setScale(scale);
        solver.setSaveBest(true);
        solver.setLowerBound(positivity ? 0.0 : Double.NEGATIVE_INFINITY);
        solver.setUpperBound(Double.POSITIVE_INFINITY);
        solver.setMaximumIterations(nbIterDeconv);
        solver.setMaximumEvaluations(2*nbIterDeconv);

        this.iterHook = iterHook;
        this.finalHook = finalHook;

    }
    public ShapedArray deconv(ShapedArray objArray){

        int iter=0;
        run = true;
        solver.setInitialSolution(objArray);

        OptimTask task = solver.start();

        while (run) {
            task = solver.getTask();
            if (task == OptimTask.ERROR) {
                System.err.format("Error: %s\n", solver.getReason());
                break;
            }
            if (task == OptimTask.NEW_X || task == OptimTask.FINAL_X) {
                if(iterHook!=null){
                    iterHook.run(this,iter++);
                }
                if (task == OptimTask.FINAL_X) {
                    break;
                }
            }
            if (task == OptimTask.WARNING) {
                break;
            }
            solver.iterate();
        }
        objArray = solver.getBestSolution().asShapedArray();
        finalHook.run(this,iter);
        return objArray;

    }
    public void abort(){
        run = false;
    }
}
