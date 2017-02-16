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
 * Compute a 3D point spread function of a wide field fluorescence microscope (WFFM)
 * <p>
 * The 3D PSF is modeled after a parameterized pupil function. It is a monochromatic
 * scalar model that defines the 3D PSF h from pupil function p.
 * Both modulus ρ(i,j) and phase φ(i,j) of pupil function p are expressed on a basis
 * of Zernike polynomials Zn.
 * <p>
 * A(z) = ρ.exp(iΦ(z)) with Φ(z) = φ + 2π( z.ψ)
 * ψ the defocus function :  ψ
 * <p>
 * <p>
 * References:
 * [1] Yves Tourneur & Eric Thiébaut, Ferreol Soulez, Loïc Denis.
 * Blind deconvolution of 3d data in wide field fluorescence microscopy.
 * <p>
 * @version
 * @author Ferréol Soulez	 <ferreol.soulez@epfl.ch>
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

    protected DoubleShapedVectorSpace defocusSpace;
    protected DoubleShapedVectorSpace phaseSpace;
    protected DoubleShapedVectorSpace modulusSpace;
    protected DoubleShapedVector modulus_coefs;  // array of Zernike coefficients that describe the modulus
    protected DoubleShapedVector phase_coefs;  // array of Zernike coefficients that describe the phase
    protected DoubleShapedVector defocus_coefs;  // array of Zernike coefficients that describe the phase


    public  abstract  void computePSF();
    public  abstract  Array3D getPSF();
    abstract protected DoubleShapedVector apply_J_modulus(ShapedVector grad);
    abstract protected DoubleShapedVector apply_J_defocus(ShapedVector grad);
    abstract protected DoubleShapedVector apply_J_phi(ShapedVector grad);

    abstract protected  void setDefocus(DoubleShapedVector defoc);
    abstract protected  void setModulus(DoubleShapedVector modulus);
    abstract protected  void setPhase(DoubleShapedVector phase);

    abstract public void freePSF();

    // abstract protected  DoubleShapedVector apply_Jacobian(DoubleShapedVector phase);

    /** Initialize the WFFM PSF model containing parameters
     *  @param psfShape shape of the PSF array
     *  @param dxy lateral pixel size
     *  @param dz axial sampling step size
     *  @param Nx number of samples along lateral X-dimension
     *  @param Ny number of samples along lateral Y-dimension
     *  @param Nz number of samples along axial Z-dimension
     *  @param radial when true use only radial zernike polynomial
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

    public DoubleShapedVector apply_Jacobian(ShapedVector grad, ShapedVectorSpace xspace){
        if(xspace ==  defocusSpace){
            return apply_J_defocus( grad);
        }else if(xspace ==  phaseSpace){
            return apply_J_phi( grad);
        }else if(xspace ==  modulusSpace){
            System.out.println("xspace modulus_coefs  "+xspace.getNumber());
            return apply_J_modulus( grad);
        }else{
            throw new IllegalArgumentException("DoubleShapedVector grad does not belong to any space");
        }
    }

    protected  void setParam(DoubleShapedVector param) {
        if(param.getOwner() ==  defocusSpace){
            setDefocus(param);
        }else if(param.getOwner() ==  phaseSpace){
            setPhase(param);
        }else if(param.getOwner() ==  modulusSpace){
            System.out.println("param modulus_coefs  "+param.getNumber());
            setModulus(param);
        }else{
            throw new IllegalArgumentException("DoubleShapedVector param does not belong to any space");
        }
    }

    public void setSingle(boolean single){
        this.single = single;
    }


}