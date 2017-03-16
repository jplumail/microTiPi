/*
 * This file is part of TiPi (a Toolkit for Inverse Problems and Imaging)
 * developed by the MitiV project.
 *
 * Copyright (c) 2014 the MiTiV project, http://mitiv.univ-lyon1.fr/
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

package microTiPi.microscopy;

import mitiv.array.Array3D;
import mitiv.base.Shape;
import mitiv.linalg.shaped.DoubleShapedVector;
import mitiv.linalg.shaped.DoubleShapedVectorSpace;
import mitiv.linalg.shaped.ShapedVector;
import mitiv.linalg.shaped.ShapedVectorSpace;

/**
 * Abstract class for to model PSF of any fluorescence microscope
 *
 * @author Ferréol
 *
 */
public abstract class MicroscopeModel
{
    protected int PState=0;   // flag to prevent useless recomputation of the PSF
    protected final static boolean NORMALIZED = true;
    protected static final double DEUXPI = 2*Math.PI;
    protected double dxy; // the lateral pixel size in meter
    protected double dz; // the axial sampling step size in meter
    protected int Nx; // number of samples along lateral X-dimension
    protected int Ny; // number of samples along lateral Y-dimension
    protected int Nz; // number of samples along axial Z-dimension
    protected boolean single = false;

    protected Shape psfShape;
    protected Array3D psf; //3D point spread function

    protected DoubleShapedVectorSpace[] parameterSpace;
    protected DoubleShapedVector[] parameterCoefs;


    /** Initialize the  PSF model containing parameters
     *  @param psfShape shape of the PSF array
     *  @param dxy lateral pixel size
     *  @param dz axial sampling step size
     *  @param single single precision flag
     */
    public MicroscopeModel(Shape psfShape,
            double dxy, double dz,
            boolean single)
    {
        this.dxy = dxy;
        this.dz = dz;

        if (psfShape.rank() !=3){
            throw new IllegalArgumentException("PSF rank should be 3");
        }
        Nx = psfShape.dimension(0);
        Ny = psfShape.dimension(1);
        Nz = psfShape.dimension(2);
        this.psfShape = psfShape;
        this.single = single;
    }


    /**
     * Launch internal routines to compute PSF
     */
    abstract public void computePSF();

    /**
     * @return the PSF
     */
    abstract public Array3D getPSF();

    /**
     * Setter for PSF parameters. The parameter type is given by the parameter space of @param
     * @param param PSF parameters
     */
    abstract public void setParam(DoubleShapedVector param);

    /**
     * Apply the Jacobian to the gradient on the PSF to get the
     *  derivative with respect to the PSF parameters
     *
     * @param grad derivative with respect to the PSF pixels
     * @param xspace PSF parameter space
     * @return derivative with respect to the PSF parameters
     */
    abstract public DoubleShapedVector apply_Jacobian(ShapedVector grad, ShapedVectorSpace xspace);


    /**
     * Free some memory
     */
    abstract public void freePSF();

    /**
     * Setter for the single precision flag
     * @param single
     */
    public void setSingle(boolean single){
        this.single = single;
    }

    /**
     * @param best_x
     * @param flag
     */
    public void set(DoubleShapedVector best_x, int flag) {
        // TODO Auto-generated method stub

    }


}