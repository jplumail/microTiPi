package microTiPi.biphoton;

import microTiPi.microscopy.MicroscopeModel;
import mitiv.array.Array3D;
import mitiv.base.Shape;
import mitiv.linalg.shaped.DoubleShapedVector;
import mitiv.linalg.shaped.ShapedVector;
import mitiv.linalg.shaped.ShapedVectorSpace;

public class BiPhotonModel extends MicroscopeModel {


    public BiPhotonModel(Shape psfShape, double NA, double lambda, double ni, double dxy, double dz, int Nx, int Ny, int Nz,
            boolean radial) {
        super(psfShape, dxy, dz, radial);
    }

    @Override
    public
    void computePSF() {
        // TODO Auto-generated method stub

    }

    /* (non-Javadoc)
     * @see microTiPi.microscopy.MicroscopeModel#getPSF()
     */
    @Override
    public Array3D getPSF() {
        // TODO Auto-generated method stub
        return null;
    }

    /* (non-Javadoc)
     * @see microTiPi.microscopy.MicroscopeModel#setParam(mitiv.linalg.shaped.DoubleShapedVector)
     */
    @Override
    public void setParam(DoubleShapedVector param) {
        // TODO Auto-generated method stub

    }

    /* (non-Javadoc)
     * @see microTiPi.microscopy.MicroscopeModel#apply_Jacobian(mitiv.linalg.shaped.ShapedVector, mitiv.linalg.shaped.ShapedVectorSpace)
     */
    @Override
    public DoubleShapedVector apply_Jacobian(ShapedVector grad, ShapedVectorSpace xspace) {
        // TODO Auto-generated method stub
        return null;
    }

    /* (non-Javadoc)
     * @see microTiPi.microscopy.MicroscopeModel#freePSF()
     */
    @Override
    public void freePSF() {
        // TODO Auto-generated method stub

    }
}

