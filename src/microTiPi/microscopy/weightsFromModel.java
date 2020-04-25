/**
 *
 */
package microTiPi.microscopy;

import microTiPi.microUtils.BlindDeconvJob;
import mitiv.array.ByteArray;
import mitiv.array.ShapedArray;
import mitiv.jobs.DeconvolutionJob;
import mitiv.utils.HistoMap;
import mitiv.utils.WeightFactory;
import mitiv.utils.WeightUpdater;

/**
 * @author ferreol
 *
 */
public class weightsFromModel implements WeightUpdater {

    private ShapedArray dataArray;
    private ByteArray badpixArray;
    private ShapedArray wghtArray;
    private double alpha, beta;

    public weightsFromModel(ShapedArray dataArray,ByteArray badpixArray) {
        this.dataArray = dataArray;
        this.badpixArray = badpixArray;
    }

    @Override
    public void update(Object caller) {
        String tmp = caller.getClass().getSimpleName();
        //        if( caller.getClass().getSimpleName() =="BlindDeconvJob") {
        if( caller  instanceof BlindDeconvJob) {
            ShapedArray modelArray = ((BlindDeconvJob) caller).getDeconvolver().getModel();
            HistoMap hm = new HistoMap(modelArray, dataArray, badpixArray);
            WeightFactory.normalize( wghtArray);
            alpha = hm.getAlpha();
            beta  = hm.getBeta();
            wghtArray= hm.computeWeightMap(modelArray);
            ((BlindDeconvJob) caller).getDeconvolver().updateWeight( wghtArray);
        }else if(  caller  instanceof DeconvolutionJob) {
            ShapedArray modelArray =  ((DeconvolutionJob) caller).getModel();
            HistoMap hm = new HistoMap(modelArray, dataArray, badpixArray);
            wghtArray= hm.computeWeightMap(modelArray);
            alpha = hm.getAlpha();
            beta  = hm.getBeta();
            WeightFactory.normalize( wghtArray);
            ((DeconvolutionJob) caller).updateWeight(wghtArray);
        }else {

        }
    }

    @Override
    public ShapedArray getWeights() {
        return wghtArray;
    }

    public double getAlpha() {
        return alpha;
    }

    public double getBeta() {
        return beta;
    }

}
