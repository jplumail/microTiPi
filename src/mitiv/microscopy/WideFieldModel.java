package mitiv.microscopy;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.jtransforms.fft.DoubleFFT_2D;
import org.jtransforms.fft.FloatFFT_2D;

import mitiv.array.Array3D;
import mitiv.array.Array4D;
import mitiv.array.Double1D;
import mitiv.array.Double2D;
import mitiv.array.Double3D;
import mitiv.array.Double4D;
import mitiv.array.Float2D;
import mitiv.array.Float3D;
import mitiv.array.Float4D;
import mitiv.base.Shape;
import mitiv.linalg.shaped.DoubleShapedVector;
import mitiv.linalg.shaped.DoubleShapedVectorSpace;
import mitiv.linalg.shaped.ShapedVector;
import mitiv.old.MathUtils;

public class WideFieldModel extends MicroscopeModel{



    protected double deltaX=0;   // position in X of the center of the defocus function inside the pupil
    protected double deltaY=0;    // position in X of the center of the defocus function inside the pupil
    protected int Nzern; // number of Zernike modes

    protected double lambda; // the emission wavelength in meters
    protected double NA; // the numerical aperture
    protected double ni; // the refractive index of the immersion medium

    protected double lambda_ni;  // (ni / \lambda)
    protected double radius; // radius of the pupil in meter^-1
    protected double pupil_area; // area of the pupil
    protected double[] Z; // Zernike polynomials basis
    protected boolean[] maskPupil; // position in the where the pupil is non null including vignetting
    protected boolean[] mapPupil; // position in the where the pupil is non null
    protected double[] rho; // pupil modulus based on Zernike polynomials
    protected double[] phi; // pupil phase based on Zernike polynomials
    protected double[] psi; // defocus function
    protected Array4D cpxPsf; // Fourier transform of the pupil function

    protected Shape cpxPsfShape;
    protected Shape aShape;
    protected Shape psf2DShape;

    // protected  Object FFT2D;

    protected int nModulus;
    protected int nDefocus;
    protected int nPhase;

    private boolean para=true;

    public WideFieldModel(Shape psfShape, double NA, double lambda, double ni, double dxy, double dz, boolean radial, boolean single){
        this( psfShape,0, 0, NA,  lambda,  ni,  dxy,  dz,  radial,  single) ;
    }

    public WideFieldModel(Shape psfShape,int nPhase, int nModulus,
            double NA, double lambda, double ni, double dxy, double dz, boolean radial, boolean single) {

        super(psfShape,    dxy, dz,  radial, single);
        if(Nx != Ny){
            throw new IllegalArgumentException("Nx should equal Ny");
        }
        this.lambda = lambda;
        this.ni = ni;
        this.Nzern = 4;
        this.NA = NA;
        this.radius = NA/lambda;
        this.lambda_ni = ni/lambda;
        this.phi = new double[Ny*Nx];
        this.psi = new double[Ny*Nx];
        cpxPsfShape = new Shape(2,Nx, Ny,  Nz);
        aShape = new Shape(2,Nx, Ny);
        psf2DShape = new Shape(Nx, Ny);

        computeMaskPupil();
        computeZernike();

        if(nModulus<1){
            nModulus = 1;
        }

        defocusSpace = new DoubleShapedVectorSpace(3);
        defocus_coefs =   defocusSpace.wrap((new double[] {ni/lambda, deltaX, deltaY}));
        setDefocus(defocus_coefs);

        if(nPhase>0){
            phaseSpace = new DoubleShapedVectorSpace(nPhase);
            phase_coefs = phaseSpace.create(0.);
            setPhase(phase_coefs);
        }

        modulusSpace =  new DoubleShapedVectorSpace(nModulus);
        modulus_coefs = modulusSpace.create(0.);
        modulus_coefs.set(0, 1.);
        setModulus(modulus_coefs);



    }


    /**
     * Compute the Zernike basis Z.
     */
    protected void computeZernike(){
        Z = Zernike.zernikeArray(Nzern, Nx, Ny, radius*dxy*Nx, NORMALIZED,radial);
        Z = MathUtils.gram_schmidt_orthonormalization(Z, Nx, Ny, Nzern);
    }

    /**
     * Compute the point spread function
     * <p>
     * h_k(z) = |a_k(z)|² = |Σ_j F_{j,k} A_j(z)|²
     */

    @Override
    public void computePSF(){


        if (PState>0)
            return;
        if(single){
            //   this.psf = new double[Nz*Ny*Nx];
            cpxPsf = Float4D.create(  cpxPsfShape);
            psf = Float3D.create( psfShape);

            final float PSFnorm = (float) (1.0/(Nx*Ny*Nz));
            final int Npix = Nx*Ny;

            int threads = Runtime.getRuntime().availableProcessors();
            ExecutorService service = Executors.newFixedThreadPool(threads);

            List<Future<GetPsfParaOut>> futures = new ArrayList<Future<GetPsfParaOut>>();
            //    ConcurrencyUtils.setNumberOfThreads(1);
            for ( int iz = 0; iz < Nz; iz++)
            {
                final int iz1 = iz;
                Callable<GetPsfParaOut> callable = new Callable<GetPsfParaOut>() {
                    @Override
                    public GetPsfParaOut call() throws Exception {
                        GetPsfParaOut output = new GetPsfParaOut(Npix,iz1,single);

                        double defoc_scale;
                        double phasePupil;
                        float[] A = new float[2*Npix];

                        if (iz1 > Nz/2)
                        {
                            defoc_scale = DEUXPI*(iz1 - Nz)*dz;
                        }
                        else
                        {
                            defoc_scale = DEUXPI*iz1*dz;
                        }

                        for (int in = 0; in < Npix; in++)
                        {
                            phasePupil = phi[in] + defoc_scale*psi[in];
                            A[2*in] = (float) (rho[in]*Math.cos(phasePupil));
                            A[2*in + 1] = (float) (rho[in]*Math.sin(phasePupil));

                        }
                        /* Fourier transform of the pupil function A(z) */
                        FloatFFT_2D FFT2D = new FloatFFT_2D(Nx, Ny);
                        FFT2D.complexForward(A);

                        for (int in = 0; in < Npix; in++)
                        {
                            ((float[])output.outA)[2*in] = A[2*in];
                            ((float[])output.outA)[2*in + 1] = -A[2*in + 1]; // store conjugate of A
                            ((float[])output.outPsf)[in] = (A[2*in]*A[2*in] + A[2*in+1]*A[2*in+1])*PSFnorm ;
                        }

                        return output;
                    }
                };
                futures.add(service.submit(callable));
            }

            service.shutdown();

            //    List<Output> outputs = new ArrayList<Output>();
            for (Future<GetPsfParaOut> future : futures) {
                GetPsfParaOut output;
                try {
                    output = future.get();
                    cpxPsf.slice(output.idxz).assign(Float3D.wrap((float[])output.outA, aShape));
                    psf.slice(output.idxz).assign(Float2D.wrap((float[])output.outPsf, psf2DShape));
                } catch (InterruptedException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                } catch (ExecutionException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                }
            }
        }else{

            if(para){
                //   this.psf = new double[Nz*Ny*Nx];
                cpxPsf = Double4D.create(  cpxPsfShape);
                psf = Double3D.create( psfShape);

                final double PSFnorm = 1.0/(Nx*Ny*Nz);
                final int Npix = Nx*Ny;

                int threads = Runtime.getRuntime().availableProcessors();
                ExecutorService service = Executors.newFixedThreadPool(threads);

                List<Future<GetPsfParaOut>> futures = new ArrayList<Future<GetPsfParaOut>>();
                //      ConcurrencyUtils.setNumberOfThreads(1);
                for ( int iz = 0; iz < Nz; iz++)
                {
                    final int iz1 = iz;
                    Callable<GetPsfParaOut> callable = new Callable<GetPsfParaOut>() {
                        @Override
                        public GetPsfParaOut call() throws Exception {
                            GetPsfParaOut output = new GetPsfParaOut(Npix,iz1,single);
                            //  ConcurrencyUtils.setNumberOfThreads(1);
                            double defoc_scale;
                            double phasePupil;
                            double[] A = new double[2*Npix];

                            if (iz1 > Nz/2)
                            {
                                defoc_scale = DEUXPI*(iz1 - Nz)*dz;
                            }
                            else
                            {
                                defoc_scale = DEUXPI*iz1*dz;
                            }

                            for (int in = 0; in < Npix; in++)
                            {
                                phasePupil = phi[in] + defoc_scale*psi[in];
                                A[2*in] = rho[in]*Math.cos(phasePupil);
                                A[2*in + 1] = rho[in]*Math.sin(phasePupil);

                            }
                            /* Fourier transform of the pupil function A(z) */


                            DoubleFFT_2D   FFT2D = new DoubleFFT_2D(Nx, Ny);

                            FFT2D.complexForward(A);

                            for (int in = 0; in < Npix; in++)
                            {
                                ((double[])output.outA)[2*in] = A[2*in];
                                ((double[])output.outA)[2*in + 1] = -A[2*in + 1]; // store conjugate of A
                                ((double[])output.outPsf)[in] = (A[2*in]*A[2*in] + A[2*in+1]*A[2*in+1])*PSFnorm ;
                            }

                            return output;
                        }
                    };
                    futures.add(service.submit(callable));
                }

                service.shutdown();

                for (Future<GetPsfParaOut> future : futures) {
                    GetPsfParaOut output;
                    try {
                        output = future.get();
                        cpxPsf.slice(output.idxz).assign(Double3D.wrap((double[])output.outA, aShape));
                        psf.slice(output.idxz).assign(Double2D.wrap((double[])output.outPsf, psf2DShape));
                    } catch (InterruptedException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    } catch (ExecutionException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    }
                }
            }else{
                cpxPsf = Double4D.create(  cpxPsfShape);
                psf = Double3D.create( psfShape);

                final double PSFnorm = 1.0/(Nx*Ny*Nz);
                final int Npix = Nx*Ny;

                DoubleFFT_2D FFT2D = new DoubleFFT_2D(Nx, Ny);

                for ( int iz = 0; iz < Nz; iz++)
                {

                    double defoc_scale;
                    double phasePupil;
                    double[] A = new double[2*Npix];

                    if (iz > Nz/2)
                    {
                        defoc_scale = DEUXPI*(iz - Nz)*dz;
                    }
                    else
                    {
                        defoc_scale = DEUXPI*iz*dz;
                    }

                    for (int in = 0; in < Npix; in++)
                    {
                        phasePupil =phi[in] + defoc_scale*psi[in];
                        A[2*in] = rho[in]*Math.cos(phasePupil);
                        A[2*in + 1] = rho[in]*Math.sin(phasePupil);

                    }
                    /* Fourier transform of the pupil function A(z) */
                    FFT2D.complexForward(A);


                    for (int iy = 0; iy < Ny; iy++){
                        for (int ix = 0; ix < Nx; ix++){
                            int in = (ix+Nx*iy);

                            ((Double4D) cpxPsf).set(0, ix, iy, iz, A[2*in]);
                            ((Double4D) cpxPsf).set(1, ix, iy, iz, -A[2*in+1]);
                            ((Double3D) psf).set(ix, iy, iz, (A[2*in]*A[2*in] + A[2*in+1]*A[2*in+1])*PSFnorm);
                        }

                    }
                }
            }



        }
        PState = 1;

    }

    /**
     * Apply the Jacobian matrix to go from  the PSF space to modulus coefficients space.
     * @param q : the gradient of some criterion in the PSF space
     * @return the gradient of this criterion in the modulus coefficients space.
     */
    @Override
    public  DoubleShapedVector apply_J_modulus( ShapedVector q)
    {
        int Ci;
        final int Npix = Nx*Ny;
        double defoc_scale = 0.;
        final double PSFNorm = 1.0/(Nx*Ny*Nz);
        //   double Aq[] = new double[2*Npix];
        final double NBeta =1./modulus_coefs.norm2();
        Double1D JRho =  Double1D.create(modulusSpace.getShape());

        JRho.fill(0.);

        if(single){
            if(para){
                int threads = Runtime.getRuntime().availableProcessors();
                ExecutorService service = Executors.newFixedThreadPool(threads);

                List<Future<ApplyJPhaOut>> futures = new ArrayList<Future<ApplyJPhaOut>>();
                //       ConcurrencyUtils.setNumberOfThreads(1);
                for ( int iz = 0; iz < Nz; iz++)
                {

                    final Float2D qz = ((Float3D) q.asShapedArray()).slice(iz);
                    final int iz1 = iz;
                    Callable<ApplyJPhaOut> callable = new Callable<ApplyJPhaOut>() {
                        @Override
                        public ApplyJPhaOut call() throws Exception {


                            double defoc_scale=0;
                            float Aq[] = new float[2*Npix];
                            ApplyJPhaOut pout = new ApplyJPhaOut( modulusSpace.getNumber());

                            if (iz1 > Nz/2)
                            {
                                defoc_scale = DEUXPI*(iz1 - Nz)*dz;
                            }
                            else
                            {
                                defoc_scale = DEUXPI*iz1*dz;
                            }
                            for (int iy = 0; iy < Ny; iy++){
                                for (int ix = 0; ix < Nx; ix++){
                                    int in = (ix+Nx*iy);
                                    float qin =  qz.get(ix, iy);
                                    Aq[2*in]=  ((Float4D) cpxPsf).get(0, ix, iy, iz1 )*qin;
                                    Aq[2*in+1]=  ((Float4D) cpxPsf).get(1, ix, iy, iz1 )*qin;
                                }

                            }

                            /* Fourier transform of the pupil function A(z) */
                            FloatFFT_2D FFT2D = new FloatFFT_2D(Nx, Ny);
                            FFT2D.complexForward(Aq);


                            for (int j = 0; j < Ny; j++)
                            {
                                for (int i = 0; i < Nx; i++)
                                {
                                    int in = i + j*Nx;
                                    if(maskPupil[in] )
                                    {
                                        double ph = phi[in] + defoc_scale*psi[in];
                                        double jin = rho[in]*(Aq[2*in]*Math.sin(ph) + Aq[2*in + 1]*Math.cos(ph));
                                        for (int k = 0; k < modulusSpace.getNumber(); k++)
                                        {
                                            int Ci= k*Npix + in;
                                            pout.grd[k] += 2*PSFNorm*jin*Z[Ci]*(1 - Math.pow(modulus_coefs.get(k)*NBeta,2))*NBeta;
                                        }
                                    }
                                }
                            }


                            return  pout;
                        }
                    };
                    futures.add(service.submit(callable));
                }

                service.shutdown();

                //    List<Output> outputs = new ArrayList<Output>();
                for (Future<ApplyJPhaOut> future : futures) {
                    ApplyJPhaOut pout;
                    try {
                        pout = future.get();
                        for (int k = 0; k < modulusSpace.getNumber(); k++)
                        {
                            JRho.set(k,JRho.get(k)+ pout.grd[k]);
                        }
                    } catch (InterruptedException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    } catch (ExecutionException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    }
                }
            }else{
                //  DoubleFFT_2D FFT2D = new DoubleFFT_2D(Ny, Nx);
                double J[] = new double[Ny*Nx];

                FloatFFT_2D FFT2D = new FloatFFT_2D(Nx, Ny);
                for (int iz = 0; iz < Nz; iz++)
                {

                    float Aq[] = new float[2*Npix];
                    if (iz > Nz/2)
                    {
                        defoc_scale = DEUXPI*(iz - Nz)*dz;
                    }
                    else
                    {
                        defoc_scale = DEUXPI*iz*dz;
                    }

                    for (int iy = 0; iy < Ny; iy++){
                        for (int ix = 0; ix < Nx; ix++){
                            int in = (ix+Nx*iy);
                            float qin =  ((Float3D) q.asShapedArray()).get(ix, iy, iz);
                            Aq[2*in]=  ((Float4D) cpxPsf).get(0, ix, iy, iz )*qin;
                            Aq[2*in+1]=  ((Float4D) cpxPsf).get(1, ix, iy, iz )*qin;
                        }

                    }

                    FFT2D.complexForward(Aq);


                    for (int in = 0; in < Npix; in++)
                    {
                        Ci = iz*Npix + in;
                        double ph = phi[in] + defoc_scale*psi[in];
                        J[in] = J[in] + Aq[2*in]*Math.cos(ph) - Aq[2*in + 1]*Math.sin(ph);
                    }

                }

                for (int k = 0; k < modulusSpace.getNumber(); k++)
                {
                    double tmp = 0;
                    for (int in = 0; in < Npix; in++)
                    {
                        Ci = k*Npix + in;
                        tmp += J[in]*Z[Ci];
                    }
                    JRho.set(k,2*PSFNorm*tmp*(1 - Math.pow(modulus_coefs.get(k)*NBeta,2))*NBeta);

                }
            }

        }else{
            if(para){
                int threads = Runtime.getRuntime().availableProcessors();
                ExecutorService service = Executors.newFixedThreadPool(threads);

                List<Future<double[]>> futures = new ArrayList<Future<double[]>>();

                for ( int iz = 0; iz < Nz; iz++)
                {

                    //    final Double2D qz = ((Double3D) q.asShapedArray()).slice(iz);
                    final double[] qz = ((Double3D) q.asShapedArray()).slice(iz).flatten();
                    final double[] Az = ((Double4D) cpxPsf).slice(iz).flatten();
                    final int iz1 = iz;
                    Callable<double[]> callable = new Callable<double[]>() {
                        @Override
                        public double[] call() throws Exception {


                            double defoc_scale=0;
                            double[] Aq = new double[2*Npix];
                            double[] J = new double[Npix];

                            if (iz1 > Nz/2)
                            {
                                defoc_scale = DEUXPI*(iz1 - Nz)*dz;
                            }
                            else
                            {
                                defoc_scale = DEUXPI*iz1*dz;
                            }
                            for (int iy = 0; iy < Ny; iy++){
                                for (int ix = 0; ix < Nx; ix++){
                                    int in = (ix+Nx*iy);
                                    double qin =  qz[in];
                                    //   Aq[2*in]=  ((Double4D) cpxPsf).get(0, ix, iy, iz1 )*qin;
                                    //  Aq[2*in+1]=  ((Double4D) cpxPsf).get(1, ix, iy, iz1 )*qin;
                                    Aq[2*in]=   Az[2*in]*qin;
                                    Aq[2*in+1]=  Az[2*in +1]*qin;
                                }

                            }

                            /* Fourier transform of the pupil function A(z) */
                            DoubleFFT_2D FFT2D = new DoubleFFT_2D(Nx, Ny);
                            FFT2D.complexForward(Aq);

                            for (int in = 0; in < Npix; in++)
                            {
                                double ph = phi[in] + defoc_scale*psi[in];
                                J[in] = J[in] + Aq[2*in]*Math.cos(ph) - Aq[2*in + 1]*Math.sin(ph);
                            }


                            /*

                            for (int j = 0; j < Ny; j++)
                            {
                                for (int i = 0; i < Nx; i++)
                                {
                                    int in = i + j*Nx;
                                    if(maskPupil[in] )
                                    {
                                        double ph = phi[in] + defoc_scale*psi[in];
                                        double jin = rho[in]*(Aq[2*in]*Math.sin(ph) + Aq[2*in + 1]*Math.cos(ph));
                                        for (int k = 0; k < modulusSpace.getNumber(); k++)
                                        {
                                            int Ci= k*Npix + in;
                                            pout.grd[k] += 2*PSFNorm*jin*Z[Ci]*(1 - Math.pow(modulus_coefs.get(k)*NBeta,2))*NBeta;

                                        }
                                    }
                                }
                            }


                            return  pout;
                             */
                            return J;
                        }
                    };
                    futures.add(service.submit(callable));
                }

                service.shutdown();
                /*
                for (Future<ApplyJPhaOut> future : futures) {
                    ApplyJPhaOut pout;
                    try {
                        pout = future.get();
                        for (int k = 0; k < modulusSpace.getNumber(); k++)
                        {
                            JRho.set(k,JRho.get(k)+ pout.grd[k]);
                        }
                    } catch (InterruptedException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    } catch (ExecutionException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    }
                }*/
                for (Future<double[]> future : futures) {
                    double[] jt;
                    try {
                        jt = future.get();
                        for (int k = 0; k < modulusSpace.getNumber(); k++)
                        {
                            double tmp = 0;
                            for (int in = 0; in < Npix; in++)
                            {
                                Ci = k*Npix + in;
                                tmp += jt[in]*Z[Ci];
                            }
                            JRho.set(k,2*PSFNorm*tmp*(1 - Math.pow(modulus_coefs.get(k)*NBeta,2))*NBeta);
                        }
                    } catch (InterruptedException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    } catch (ExecutionException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    }
                }
            }else{
                //     DoubleFFT_2D FFT2D = new DoubleFFT_2D(Ny, Nx);
                double J[] = new double[Ny*Nx];


                double Aq[] = new double[2*Npix];
                DoubleFFT_2D FFT2D = new DoubleFFT_2D(Nx, Ny);
                for (int iz = 0; iz < Nz; iz++)
                {

                    if (iz > Nz/2)
                    {
                        defoc_scale = DEUXPI*(iz - Nz)*dz;
                    }
                    else
                    {
                        defoc_scale = DEUXPI*iz*dz;
                    }

                    for (int iy = 0; iy < Ny; iy++){
                        for (int ix = 0; ix < Nx; ix++){
                            int in = (ix+Nx*iy);
                            double qin =  ((Double3D) q.asShapedArray()).get(ix, iy, iz);
                            Aq[2*in]=  ((Double4D) cpxPsf).get(0, ix, iy, iz )*qin;
                            Aq[2*in+1]=  ((Double4D) cpxPsf).get(1, ix, iy, iz )*qin;
                        }

                    }


                    FFT2D.complexForward(Aq);

                    for (int in = 0; in < Npix; in++)
                    {
                        Ci = iz*Npix + in;
                        double ph = phi[in] + defoc_scale*psi[in];
                        J[in] = J[in] + Aq[2*in]*Math.cos(ph) - Aq[2*in + 1]*Math.sin(ph);
                    }

                }


                for (int k = 0; k < modulusSpace.getNumber(); k++)
                {
                    double tmp = 0;
                    for (int in = 0; in < Npix; in++)
                    {
                        Ci = k*Npix + in;
                        tmp += J[in]*Z[Ci];
                    }
                    JRho.set(k,2*PSFNorm*tmp*(1 - Math.pow(modulus_coefs.get(k)*NBeta,2))*NBeta);
                    //  JRho.set(k,2*PSFNorm*tmp);
                }
            }
        }

        return modulusSpace.create(JRho);

    }


    /**
     * Apply the Jacobian matrix to go from  the PSF space to phase coefficients space.
     * @param q : the gradient of some criterion in the PSF space
     * @return the gradient of this criterion in the phase coefficients space.
     */
    @Override
    public DoubleShapedVector apply_J_phi(ShapedVector q)
    {
        int Ci;
        final int Npix = Nx*Ny;
        final double PSFNorm = 1.0/(Nx*Ny*Nz);
        Double1D JPhi =  Double1D.create(phaseSpace.getShape());
        JPhi.fill(0.);
        if(single){
            if(para){
                int threads = Runtime.getRuntime().availableProcessors();
                ExecutorService service = Executors.newFixedThreadPool(threads);

                List<Future<ApplyJPhaOut>> futures = new ArrayList<Future<ApplyJPhaOut>>();
                //        ConcurrencyUtils.setNumberOfThreads(1);
                for ( int iz = 0; iz < Nz; iz++)
                {

                    final Float2D qz = ((Float3D) q.asShapedArray()).slice(iz);
                    final int iz1 = iz;
                    Callable<ApplyJPhaOut> callable = new Callable<ApplyJPhaOut>() {
                        @Override
                        public ApplyJPhaOut call() throws Exception {


                            double defoc_scale=0;
                            float Aq[] = new float[2*Npix];
                            ApplyJPhaOut pout = new ApplyJPhaOut( phaseSpace.getNumber());

                            if (iz1 > Nz/2)
                            {
                                defoc_scale = DEUXPI*(iz1 - Nz)*dz;
                            }
                            else
                            {
                                defoc_scale = DEUXPI*iz1*dz;
                            }
                            for (int iy = 0; iy < Ny; iy++){
                                for (int ix = 0; ix < Nx; ix++){
                                    int in = (ix+Nx*iy);
                                    float qin =  qz.get(ix, iy);
                                    Aq[2*in]=  ((Float4D) cpxPsf).get(0, ix, iy, iz1 )*qin;
                                    Aq[2*in+1]=  ((Float4D) cpxPsf).get(1, ix, iy, iz1 )*qin;
                                }

                            }

                            /* Fourier transform of the pupil function A(z) */
                            FloatFFT_2D FFT2D = new FloatFFT_2D(Nx, Ny);
                            FFT2D.complexForward(Aq);


                            for (int j = 0; j < Ny; j++)
                            {
                                for (int i = 0; i < Nx; i++)
                                {
                                    int in = i + j*Nx;
                                    if(maskPupil[in] )
                                    {
                                        double ph = phi[in] + defoc_scale*psi[in];
                                        double jin = rho[in]*(Aq[2*in]*Math.sin(ph) + Aq[2*in + 1]*Math.cos(ph));
                                        for (int k = 0; k < phaseSpace.getNumber(); k++)
                                        {
                                            int Ci;

                                            if(radial){
                                                Ci= (k+1)*Npix + in;
                                            }else{
                                                Ci= (k+3)*Npix + in;
                                            }
                                            pout.grd[k] -= 2*PSFNorm*jin*Z[Ci];
                                        }
                                    }
                                }
                            }


                            return  pout;
                        }
                    };
                    futures.add(service.submit(callable));
                }

                service.shutdown();

                for (Future<ApplyJPhaOut> future : futures) {
                    ApplyJPhaOut pout;
                    try {
                        pout = future.get();
                        for (int k = 0; k < phaseSpace.getNumber(); k++)
                        {
                            JPhi.set(k,JPhi.get(k)+ pout.grd[k]);
                        }
                    } catch (InterruptedException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    } catch (ExecutionException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    }
                }
            }else{
                double J[] = new double[Ny*Nx];
                float[] Aq = new float[2*Npix];
                FloatFFT_2D FFT2D = new FloatFFT_2D(Nx, Ny);
                for (int iz = 0; iz < Nz; iz++)
                {

                    double defoc_scale=0.;
                    if (iz > Nz/2)
                    {
                        defoc_scale = DEUXPI*(iz - Nz)*dz;
                    }
                    else
                    {
                        defoc_scale = DEUXPI*iz*dz;
                    }

                    for (int iy = 0; iy < Ny; iy++){
                        for (int ix = 0; ix < Nx; ix++){
                            int in = (ix+Nx*iy);
                            float qin =  ((Float3D) q.asShapedArray()).get(ix, iy, iz);
                            Aq[2*in]=  ((Float4D) cpxPsf).get(0, ix, iy, iz )*qin;
                            Aq[2*in+1]=  ((Float4D) cpxPsf).get(1, ix, iy, iz )*qin;
                        }

                    }

                    FFT2D.complexForward(Aq);

                    for (int in = 0; in < Npix; in++)
                    {
                        Ci = iz*Npix + in;
                        double ph = phi[in] + defoc_scale*psi[in];
                        J[in] = J[in] + rho[in]*(Aq[2*in]*Math.sin(ph) + Aq[2*in + 1]*Math.cos(ph));
                    }
                }


                for (int k = 0; k < phaseSpace.getNumber(); k++)
                {
                    double tmp = 0;
                    for (int in = 0; in < Npix; in++)
                    {
                        Ci = k*Npix + in;
                        if(radial){
                            tmp += J[in]*Z[Ci + 1*Npix];
                        }else{
                            tmp += J[in]*Z[Ci + 3*Npix];
                        }
                    }
                    JPhi.set(k, -2*PSFNorm*tmp);
                }
            }
        }else{
            if(para){
                int threads = Runtime.getRuntime().availableProcessors();
                ExecutorService service = Executors.newFixedThreadPool(threads);

                List<Future<ApplyJPhaOut>> futures = new ArrayList<Future<ApplyJPhaOut>>();
                //        ConcurrencyUtils.setNumberOfThreads(1);
                for ( int iz = 0; iz < Nz; iz++)
                {

                    final Double2D qz = ((Double3D) q.asShapedArray()).slice(iz);
                    final int iz1 = iz;
                    Callable<ApplyJPhaOut> callable = new Callable<ApplyJPhaOut>() {
                        @Override
                        public ApplyJPhaOut call() throws Exception {


                            double defoc_scale=0;
                            double Aq[] = new double[2*Npix];
                            ApplyJPhaOut pout = new ApplyJPhaOut( phaseSpace.getNumber());

                            if (iz1 > Nz/2)
                            {
                                defoc_scale = DEUXPI*(iz1 - Nz)*dz;
                            }
                            else
                            {
                                defoc_scale = DEUXPI*iz1*dz;
                            }
                            for (int iy = 0; iy < Ny; iy++){
                                for (int ix = 0; ix < Nx; ix++){
                                    int in = (ix+Nx*iy);
                                    double qin =  qz.get(ix, iy);
                                    Aq[2*in]=  ((Double4D) cpxPsf).get(0, ix, iy, iz1 )*qin;
                                    Aq[2*in+1]=  ((Double4D) cpxPsf).get(1, ix, iy, iz1 )*qin;
                                }

                            }

                            /* Fourier transform of the pupil function A(z) */
                            DoubleFFT_2D FFT2D = new DoubleFFT_2D(Nx, Ny);
                            FFT2D.complexForward(Aq);

                            for (int j = 0; j < Ny; j++)
                            {
                                for (int i = 0; i < Nx; i++)
                                {
                                    int in = i + j*Nx;
                                    if(maskPupil[in] )
                                    {
                                        double ph = phi[in] + defoc_scale*psi[in];
                                        double jin = rho[in]*(Aq[2*in]*Math.sin(ph) + Aq[2*in + 1]*Math.cos(ph));
                                        for (int k = 0; k < phaseSpace.getNumber(); k++)
                                        {
                                            int Ci;

                                            if(radial){
                                                Ci= (k+1)*Npix + in;
                                            }else{
                                                Ci= (k+3)*Npix + in;
                                            }
                                            pout.grd[k] -= 2*PSFNorm*jin*Z[Ci];
                                        }
                                    }
                                }
                            }


                            return  pout;
                        }
                    };
                    futures.add(service.submit(callable));
                }

                service.shutdown();

                for (Future<ApplyJPhaOut> future : futures) {
                    ApplyJPhaOut pout;
                    try {
                        pout = future.get();
                        for (int k = 0; k < phaseSpace.getNumber(); k++)
                        {
                            JPhi.set(k,JPhi.get(k)+ pout.grd[k]);
                        }
                    } catch (InterruptedException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    } catch (ExecutionException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    }
                }
            }else{
                double J[] = new double[Ny*Nx];
                double[] Aq = new double[2*Npix];
                DoubleFFT_2D FFT2D = new DoubleFFT_2D(Ny, Nx);
                for (int iz = 0; iz < Nz; iz++)
                {

                    double defoc_scale=0.;
                    if (iz > Nz/2)
                    {
                        defoc_scale = DEUXPI*(iz - Nz)*dz;
                    }
                    else
                    {
                        defoc_scale = DEUXPI*iz*dz;
                    }

                    for (int iy = 0; iy < Ny; iy++){
                        for (int ix = 0; ix < Nx; ix++){
                            int in = (ix+Nx*iy);
                            double qin =  ((Double3D) q.asShapedArray()).get(ix, iy, iz);
                            Aq[2*in]=  ((Double4D) cpxPsf).get(0, ix, iy, iz )*qin;
                            Aq[2*in+1]=  ((Double4D) cpxPsf).get(1, ix, iy, iz )*qin;
                        }

                    }

                    FFT2D.complexForward(Aq);

                    for (int in = 0; in < Npix; in++)
                    {
                        Ci = iz*Npix + in;
                        double ph = phi[in] + defoc_scale*psi[in];
                        J[in] = J[in] + rho[in]*(Aq[2*in]*Math.sin(ph) + Aq[2*in + 1]*Math.cos(ph));
                    }
                }


                for (int k = 0; k < phaseSpace.getNumber(); k++)
                {
                    double tmp = 0;
                    for (int in = 0; in < Npix; in++)
                    {
                        Ci = k*Npix + in;
                        if(radial){
                            tmp += J[in]*Z[Ci + 1*Npix];
                        }else{
                            tmp += J[in]*Z[Ci + 3*Npix];
                        }
                    }
                    JPhi.set(k, -2*PSFNorm*tmp);
                }
            }
        }
        return phaseSpace.create(JPhi );
    }


    /**
     * Apply the Jacobian matrix to go from  the PSF space to defocus coefficients space.
     * @param q : the gradient of some criterion in the PSF space
     * @return the gradient of this criterion in the defocus coefficients space.
     */
    @Override
    public DoubleShapedVector apply_J_defocus(ShapedVector q)
    {

        // long startTime = System.currentTimeMillis();

        double scale_x = 1/(Nx*dxy);
        double scale_y = 1/(Ny*dxy);
        //  double defoc, tmpvar, idef;
        double d0 = 0, d1 = 0, d2 = 0;
        final double[] rx = new double[Nx];
        final double[] ry = new double[Ny];
        final int Npix =  Nx*Ny;
        final double PSFNorm = 1.0/(Nx*Ny*Nz);
        double[] grd = new double[defocusSpace.getNumber()];

        for(int nx = 0; nx < Nx; nx++)
        {
            if(nx > Nx/2)
            {
                rx[nx] = (nx - Nx)*scale_x - deltaX;
            }
            else
            {
                rx[nx]  = nx*scale_x - deltaX;
            }
        }

        for(int ny = 0; ny < Ny; ny++)
        {
            if(ny > Ny/2)
            {
                ry[ny] = (ny - Ny)*scale_y - deltaY;
            }else
            {
                ry[ny]  = ny*scale_y - deltaY;
            }
        }
        if(single){
            if(para){
                int threads = Runtime.getRuntime().availableProcessors();
                ExecutorService service = Executors.newFixedThreadPool(threads);

                List<Future<ApplyJDefOut>> futures = new ArrayList<Future<ApplyJDefOut>>();
                //  ConcurrencyUtils.setNumberOfThreads(1);
                for ( int iz = 0; iz < Nz; iz++)
                {

                    final Float2D qz = ((Float3D) q.asShapedArray()).slice(iz);
                    final int iz1 = iz;
                    Callable<ApplyJDefOut> callable = new Callable<ApplyJDefOut>() {
                        @Override
                        public ApplyJDefOut call() throws Exception {


                            double defoc_scale=0;
                            float Aq[] = new float[2*Npix];
                            double defoc;
                            ApplyJDefOut dout = new ApplyJDefOut(0,0,0);
                            if (iz1 > Nz/2)
                            {
                                defoc = (iz1 - Nz)*dz;
                                defoc_scale = DEUXPI*(iz1 - Nz)*dz;
                            }
                            else
                            {
                                defoc = iz1*dz;
                                defoc_scale = DEUXPI*iz1*dz;
                            }

                            for (int iy = 0; iy < Ny; iy++){
                                for (int ix = 0; ix < Nx; ix++){
                                    int in = (ix+Nx*iy);
                                    float qin =  qz.get(ix, iy);
                                    Aq[2*in]=  ((Float4D) cpxPsf).get(0, ix, iy, iz1 )*qin;
                                    Aq[2*in+1]=  ((Float4D) cpxPsf).get(1, ix, iy, iz1 )*qin;
                                }

                            }
                            /* Fourier transform of the pupil function A(z) */
                            FloatFFT_2D FFT2D = new FloatFFT_2D(Nx, Ny);
                            FFT2D.complexForward(Aq);

                            for (int j = 0; j < Ny; j++)
                            {
                                for (int i = 0; i < Nx; i++)
                                {
                                    int in = i + j*Nx;
                                    if(maskPupil[in] )
                                    {
                                        double idef= 1./psi[in];
                                        double ph = phi[in] + defoc_scale*psi[in];
                                        double tmpvar = -DEUXPI*rho[in]*( Aq[2*in]*Math.sin(ph) + Aq[2*in + 1]*Math.cos(ph) )*PSFNorm;
                                        {
                                            dout.d1 -= tmpvar*( rx[i]*(defoc*idef ));
                                            dout.d2 -= tmpvar*( ry[j]*(defoc*idef) );
                                            dout.d0 += tmpvar*( idef*lambda_ni*defoc );
                                        }
                                    }
                                }
                            }

                            return  dout;
                        }
                    };
                    futures.add(service.submit(callable));
                }

                service.shutdown();

                for (Future<ApplyJDefOut> future : futures) {
                    ApplyJDefOut output;
                    try {
                        output = future.get();
                        d0 += output.d0;
                        d1 -= output.d1;
                        d2 -= output.d2;
                    } catch (InterruptedException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    } catch (ExecutionException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    }
                }
            }else{

                FloatFFT_2D FFT2D = new FloatFFT_2D(Nx, Ny);
                double defoc, idef, tmpvar;
                float Aq[] = new float[2*Npix];

                for (int iz = 0; iz < Nz; iz++)
                {
                    double defoc_scale =0.;
                    if (iz > Nz/2)
                    {
                        defoc = (iz - Nz)*dz;
                        defoc_scale  = DEUXPI*defoc;
                    }
                    else
                    {
                        defoc = iz*dz;
                        defoc_scale = DEUXPI*defoc;
                    }


                    for (int iy = 0; iy < Ny; iy++){
                        for (int ix = 0; ix < Nx; ix++){
                            int in = (ix+Nx*iy);
                            float qin =  ((Float3D) q.asShapedArray()).get(ix, iy, iz);
                            Aq[2*in]=  ((Float4D) cpxPsf).get(0, ix, iy, iz )*qin;
                            Aq[2*in+1]=  ((Float4D) cpxPsf).get(1, ix, iy, iz )*qin;
                        }

                    }


                    FFT2D.complexForward(Aq);



                    for (int j = 0; j < Ny; j++)
                    {
                        for (int i = 0; i < Nx; i++)
                        {
                            int in = i + j*Nx;
                            if(maskPupil[in] )
                            {
                                idef= 1./psi[in];
                                double ph = phi[in] + defoc_scale*psi[in];
                                tmpvar = -DEUXPI*rho[in]*( Aq[2*in]*Math.sin(ph) + Aq[2*in + 1]*Math.cos(ph) )*PSFNorm;
                                {
                                    d1 -= tmpvar*( rx[i]*(defoc*idef ));
                                    d2 -= tmpvar*( ry[j]*(defoc*idef) );
                                    d0 += tmpvar*( idef*lambda_ni*defoc );
                                }
                            }
                        }
                    }
                }

            }
        }else{
            if(para){
                int threads = Runtime.getRuntime().availableProcessors();
                ExecutorService service = Executors.newFixedThreadPool(threads);

                List<Future<double[]>> futures = new ArrayList<Future<double[]>>();
                //  ConcurrencyUtils.setNumberOfThreads(1);
                for ( int iz = 0; iz < Nz; iz++)
                {
                    final double[] qz = ((Double3D) q.asShapedArray()).slice(iz).flatten();
                    final double[] Az = ((Double4D) cpxPsf).slice(iz).flatten();

                    //      final Double2D qz = ((Double3D) q.asShapedArray()).slice(iz);
                    final int iz1 = iz;
                    Callable<double[]> callable = new Callable<double[]>() {
                        @Override
                        public double[] call() throws Exception {


                            double defoc_scale=0;
                            double Aq[] = new double[2*Npix];
                            double defoc;
                            //  ApplyJDefOut dout = new ApplyJDefOut(0,0,0);
                            double[] dout = new double[3];
                            if (iz1 > Nz/2)
                            {
                                defoc = (iz1 - Nz)*dz;
                                defoc_scale = DEUXPI*(iz1 - Nz)*dz;
                            }
                            else
                            {
                                defoc = iz1*dz;
                                defoc_scale = DEUXPI*iz1*dz;
                            }

                            for (int iy = 0; iy < Ny; iy++){
                                for (int ix = 0; ix < Nx; ix++){
                                    int in = (ix+Nx*iy);
                                    //  double qin =  qz.get(ix, iy);
                                    //  Aq[2*in]=  ((Double4D) cpxPsf).get(0, ix, iy, iz1 )*qin;
                                    //  Aq[2*in+1]=  ((Double4D) cpxPsf).get(1, ix, iy, iz1 )*qin;
                                    Aq[2*in]=   Az[2*in]*qz[in];
                                    Aq[2*in+1]=  Az[2*in+1]*qz[in];
                                }

                            }
                            /* Fourier transform of the pupil function A(z) */
                            DoubleFFT_2D FFT2D = new DoubleFFT_2D(Nx, Ny);
                            FFT2D.complexForward(Aq);

                            for (int j = 0; j < Ny; j++)
                            {
                                for (int i = 0; i < Nx; i++)
                                {
                                    int in = i + j*Nx;
                                    if(maskPupil[in] )
                                    {
                                        double idef= 1./psi[in];
                                        double ph = phi[in] + defoc_scale*psi[in];
                                        double tmpvar = -DEUXPI*rho[in]*( Aq[2*in]*Math.sin(ph) + Aq[2*in + 1]*Math.cos(ph) )*PSFNorm;
                                        {
                                            /* dout.d1 -= tmpvar*( rx[i]*(defoc*idef ));
                                            dout.d2 -= tmpvar*( ry[j]*(defoc*idef) );
                                            dout.d0 += tmpvar*( idef*lambda_ni*defoc );*/
                                            dout[1] -= tmpvar*( rx[i]*(defoc*idef ));
                                            dout[2] -= tmpvar*( ry[j]*(defoc*idef) );
                                            dout[0] += tmpvar*( idef*lambda_ni*defoc );
                                        }
                                    }
                                }
                            }

                            return  dout;
                        }
                    };
                    futures.add(service.submit(callable));
                }

                service.shutdown();

                for (Future<double[]> future : futures) {
                    double[] output;
                    try {
                        output = future.get();
                        d0 += output[0];
                        d1 -= output[1];
                        d2 -= output[2];
                    } catch (InterruptedException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    } catch (ExecutionException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    }
                }
            }else{

                DoubleFFT_2D FFT2D = new DoubleFFT_2D(Nx, Ny);
                for ( int iz = 0; iz < Nz; iz++)
                {
                    final double[] qz = ((Double3D) q.asShapedArray()).slice(iz).flatten();
                    final double[] Az = ((Double4D) cpxPsf).slice(iz).flatten();

                    //      final Double2D qz = ((Double3D) q.asShapedArray()).slice(iz);
                    final int iz1 = iz;


                    double defoc_scale=0;
                    double Aq[] = new double[2*Npix];
                    double defoc;
                    //  ApplyJDefOut dout = new ApplyJDefOut(0,0,0);
                    // double[] dout = new double[3];
                    if (iz1 > Nz/2)
                    {
                        defoc = (iz1 - Nz)*dz;
                        defoc_scale = DEUXPI*(iz1 - Nz)*dz;
                    }
                    else
                    {
                        defoc = iz1*dz;
                        defoc_scale = DEUXPI*iz1*dz;
                    }

                    for (int iy = 0; iy < Ny; iy++){
                        for (int ix = 0; ix < Nx; ix++){
                            int in = (ix+Nx*iy);
                            //  double qin =  qz.get(ix, iy);
                            //  Aq[2*in]=  ((Double4D) cpxPsf).get(0, ix, iy, iz1 )*qin;
                            //  Aq[2*in+1]=  ((Double4D) cpxPsf).get(1, ix, iy, iz1 )*qin;
                            Aq[2*in]=   Az[2*in]*qz[in];
                            Aq[2*in+1]=  Az[2*in+1]*qz[in];
                        }

                    }
                    /* Fourier transform of the pupil function A(z) */
                    FFT2D.complexForward(Aq);


                    for (int j = 0; j < Ny; j++)
                    {
                        for (int i = 0; i < Nx; i++)
                        {
                            int in = i + j*Nx;
                            if(maskPupil[in] )
                            {
                                double idef= 1./psi[in];
                                double ph = phi[in] + defoc_scale*psi[in];
                                double tmpvar = -DEUXPI*rho[in]*( Aq[2*in]*Math.sin(ph) + Aq[2*in + 1]*Math.cos(ph) )*PSFNorm;
                                {
                                    /* dout.d1 -= tmpvar*( rx[i]*(defoc*idef ));
                                            dout.d2 -= tmpvar*( ry[j]*(defoc*idef) );
                                            dout.d0 += tmpvar*( idef*lambda_ni*defoc );*/
                                    d1 -= tmpvar*( rx[i]*(defoc*idef ));
                                    d2 -= tmpvar*( ry[j]*(defoc*idef) );
                                    d0 += tmpvar*( idef*lambda_ni*defoc );
                                }
                            }
                        }
                    }

                }
            }
        }
        switch(defocusSpace.getNumber())
        {
            case 3:
                grd[2] = d2;
                grd[1] = d1;
            case 1:
                grd[0] = d0;
                break;
            case 2:
                grd[2] = d2;
                grd[1] = d1;
                break;
        }


        return  defocusSpace.create(Double1D.wrap(grd, defocusSpace.getShape()));

    }


    /** Determine the map where the pupil in non null.  It sets  maskPupil and its area pupil_area
     */
    private void computeMaskPupil()
    {
        maskPupil = new boolean[Nx*Ny];
        mapPupil = new boolean[Nx*Ny];
        double scale_y = Math.pow(1/dxy/Ny, 2);
        double scale_x = Math.pow(1/dxy/Nx, 2);
        double rx, ry, ix, iy;
        double radius2 = radius*radius;
        pupil_area =0.;
        for(int ny = 0; ny < Ny; ny++)
        {
            iy = Math.min(ny, Ny - ny);
            ry = iy*iy*scale_y;
            for(int nx = 0; nx < Nx; nx++)
            {
                ix = Math.min(nx, Nx - nx);
                rx = ix*ix*scale_x;
                if( (rx + ry) < radius2 )
                {
                    maskPupil[nx + ny*Nx] = true;
                    mapPupil[nx + ny*Nx] = true;
                    pupil_area += 1;

                }else{

                    maskPupil[nx + ny*Nx] = false;
                    mapPupil[nx + ny*Nx] = false;
                }
            }
        }
        pupil_area = Math.sqrt(pupil_area);
        freePSF();
    }


    private class ApplyJDefOut {
        double d0;
        double d1;
        double d2;
        public ApplyJDefOut(double d0_,double d1_, double d2_){
            d0 = d0_;
            d1 = d1_;
            d2 = d2_;
        }
    }


    private class ApplyJPhaOut {
        double[] grd;
        public ApplyJPhaOut( int nMode){
            grd = new double[nMode];
            Arrays.fill(grd, 0.);
        }
    }


    private class GetPsfParaOut{
        Object outA;
        Object  outPsf;
        int idxz;
        public GetPsfParaOut(int nPix, int iz, boolean single){
            idxz = iz;
            if(single){
                outA = new float[2*nPix];
                outPsf = new float[2*nPix];
            }else{
                outA = new double[2*nPix];
                outPsf = new double[2*nPix];

            }
        }
    }

    /**
     * Compute the defocus aberration ψ of the phase pupil
     * <p>
     */

    public void computeDefocus()
    {
        double lambda_ni2 = lambda_ni*lambda_ni;
        double scale_x = 1/(Nx*dxy);
        double scale_y = 1/(Ny*dxy);
        double q, rx, ry;
        for (int ny = 0; ny < Ny; ny++)
        {
            if(ny > Ny/2)
            {
                ry = Math.pow(scale_y*(ny - Ny) - deltaY, 2);
            }
            else
            {
                ry = Math.pow(scale_y*ny - deltaY, 2);
            }

            for (int nx = 0; nx < Nx; nx++)
            {
                int nxy = nx + ny*Nx;
                if (mapPupil[nxy] )
                {
                    if(nx > Nx/2)
                    {
                        rx = Math.pow(scale_x*(nx - Nx) - deltaX, 2);
                    }
                    else
                    {
                        rx = Math.pow(scale_x*nx - deltaX, 2);
                    }

                    q = lambda_ni2 - rx - ry;

                    if (q < 0.0)
                    {
                        psi[nxy] = 0;
                        maskPupil[nxy] = false;
                    }
                    else
                    {
                        psi[nxy] = Math.sqrt(q);
                        maskPupil[nxy] = true;
                    }
                }
            }
        }
        freePSF();
    }



    /**
     * @param defocus Update the defocus and the depth functions according the parameters
     * defocus. Depending on the number of elements of defocus:
     * 3 :  defocus = {n_i / \lambda, \delta_x, \delta_y}
     * 2 :  defocus = { \delta_x, \delta_y}
     * 1 :  defocus = {n_i / \lambda}
     */
    @Override
    public void setDefocus(DoubleShapedVector defoc) {
        if(defoc.belongsTo(defocusSpace)){
            defocus_coefs = defoc;
        }else{
            throw new IllegalArgumentException("defocus  does not belong to the defocusSpace");
        }
        switch (defoc.getNumber())
        {
            case 3:
                deltaX = defoc.get(1);
                deltaY = defoc.get(2);
            case 1:
                lambda_ni = defoc.get(0);
                break;
            case 2:
                deltaX = defoc.get(1);
                deltaY = defoc.get(2);
                break;
            default:
                throw new IllegalArgumentException("bad defocus  parameters");
        }
        computeDefocus();
        freePSF();
    }
    /**
     * Compute the modulus ρ on a Zernike polynomial basis
     * <p>
     * The coefficients β are normalized and the modulus is
     * ρ = Σ_n β_n Z_n
     * @param beta Zernike coefficients
     */
    @Override
    public void setModulus(DoubleShapedVector beta) {
        if(beta.belongsTo(modulusSpace)){
            modulus_coefs = beta;
        }else{
            throw new IllegalArgumentException("DoubleShapedVector beta does not belong to the modulus space");
        }

        int Npix = Nx*Ny;
        rho = new double[Npix];
        //    double betaNorm = 1./(Math.sqrt(MathUtils.innerProd(modulus_coefs, modulus_coefs)));
        double betaNorm = 1./( beta.norm2());
        for(int in = 0; in < Npix; in++)
        {
            if (maskPupil[in]  )
            {

                for (int n = 0; n < beta.getNumber(); n++)
                {
                    rho[in] += Z[in + n*Npix]*beta.get(n)*betaNorm;
                }


            }
        }


        freePSF();
    }

    @Override
    public void setPhase(DoubleShapedVector phase) {
        if(phase.belongsTo(phaseSpace)){
            phase_coefs = phase;
        }else{
            throw new IllegalArgumentException("phase parameter does not belong to the right space  ");
        }
        int Npix = Nx*Ny;
        phi = new double[Npix];
        for(int in = 0; in < Npix; in++)
        {
            if (maskPupil[in] )
            {

                for (int n = 0; n < phase.getNumber(); ++n)
                {
                    if(radial){
                        phi[in] += Z[in + (n + 1)*Npix]*phase.get(n);
                    }else{
                        phi[in] += Z[in + (n + 3)*Npix]*phase.get(n);
                    }
                }
            }
        }


        freePSF();

    }






    /**
     * @return the modulus of the pupil
     */
    public double[] getRho() {
        if (PState<1){
            computePSF();
        }
        return rho;
    }


    /**
     * @return the wavelength used in the computation
     */
    public double getLambda() {
        return lambda;
    }

    /**
     * @return the refractive index of the immersion medium used in the computation
     */
    public double getNi() {
        return ni;
    }


    /**
     * @return the phase of the pupil
     */
    public double[] getPhi(){
        if (PState<1){
            computePSF();
        }
        return phi;
    }

    /**
     * @return the defocus function
     */
    public double[] getPsi() {
        if (PState<1){
            computePSF();
        }
        return psi;
    }




    /**
     * @return modulus coefficients
     */
    public DoubleShapedVector getBeta() {
        return modulus_coefs;
    }

    /**
     * @return phase coefficients
     */
    public DoubleShapedVector getAlpha() {
        return phase_coefs;
    }

    /**
     * @return defocus coefficients in 1./wavelength
     */
    public double[] getDefocusMultiplyByLambda() {
        if (PState<1){
            computePSF();
        }
        double[] defocus = {lambda_ni*lambda, deltaX*lambda, deltaY*lambda};
        return defocus;
    }

    /**
     * @return defocus coefficients
     */
    public double[] getDefocus() {
        if (PState<1){
            computePSF();
        }
        double[] defocus = {lambda_ni, deltaX, deltaY};
        return defocus;
    }


    /**
     * @return the pupil mask
     */
    public boolean[] getMaskPupil() {
        if (PState<1){
            computePSF();
        }
        return maskPupil;
    }




    /**
     * @return the PSF
     */
    @Override
    public  Array3D getPSF() {

        if (PState<1){
            computePSF();
        }
        return psf;
    }


    /**
     * @return the Zernike basis
     */
    public double[] getZernike() {
        return Z;
    }

    /**
     * @return the number of zernike polynomial used in the Zernike basis
     */
    public int getNZern() {
        return Nzern;
    }

    /**
     * @param k
     * @return the k-th zernike of the basis
     */
    public double[] getZernike(int k) {
        return MathUtils.getArray(Z, Nx, Ny, k);
    }

    /**
     * @return the complex PSF
     */
    public Array4D get_cpxPsf() {
        if (PState<1){
            computePSF();
        }
        return cpxPsf;
    }

    /**
     * Plot some information about the WideFieldModel object for debugging purpose
     */
    public void getInfo()
    {
        System.out.println("----PSF----");
        MathUtils.stat(psf.toDouble().getData());
        System.out.println();

        System.out.println("----PHI----");
        MathUtils.stat(phi);
        System.out.println();

        System.out.println("----RHO----");
        MathUtils.stat(rho);
        System.out.println();

        System.out.println("----PSI----");
        MathUtils.stat(psi);
        System.out.println();

        /*System.out.println("----PUPIL----");
        MathUtils.stat( maskPupil);
        System.out.println();*/

        System.out.println("----a----");
        MathUtils.statC(cpxPsf.toDouble().getData());
        System.out.println();

        System.out.println("----ZERNIKES----");
        MathUtils.stat(Z);
    }

    public void setNPhase(int nPhase) {
        if(nPhase>0){
            phaseSpace = new DoubleShapedVectorSpace(nPhase);
            if(radial){
                Nzern = Math.max(nPhase+1, modulusSpace.getNumber());
            }else{
                Nzern = Math.max(nPhase+3, modulusSpace.getNumber());
            }
            computeZernike();
            phase_coefs = phaseSpace.create(0.);
            setPhase(phase_coefs);
        }else{
            phaseSpace = null;
        }
    }
    public void setNModulus(int nModulus) {
        if(nModulus<1){
            nModulus = 1;
        }

        modulusSpace =  new DoubleShapedVectorSpace(nModulus);

        if(radial){
            Nzern = Math.max(phaseSpace.getNumber()+1, nModulus);
        }else{
            Nzern = Math.max(phaseSpace.getNumber()+3, nModulus);
        }
        computeZernike();
        modulus_coefs = modulusSpace.create(0.);
        modulus_coefs.set(0, 1.);

        setModulus(modulus_coefs);


    }



    /**
     * reset PSF and its complex a  to free some memory
     * Set the flag PState to 0
     */
    @Override
    public void freePSF() {
        PState =0;
        cpxPsf = null;
        psf = null;
    }



    public int getNModulus() {
        return modulus_coefs.getNumber();
    }

    public int getNPhase() {
        return phase_coefs.getNumber();
    }
}

