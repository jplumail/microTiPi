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

import mitiv.array.Array1D;
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
import mitiv.old.MathUtils;
import pl.edu.icm.jlargearrays.ConcurrencyUtils;

public class WideFieldModel extends MicroscopeModel{



    protected double deltaX=0;   // position in X of the center of the defocus function inside the pupil
    protected double deltaY=0;    // position in X of the center of the defocus function inside the pupil
    protected int Nzern; // number of Zernike modes
    protected double lambda; // the emission wavelength in meters

    protected double lambda_ni;  // (ni / \lambda)
    protected double radius; // radius of the pupil in meter
    protected double pupil_area; // area of the pupil
    protected double[] Z; // Zernike polynomials basis
    protected double[] maskPupil; // position in the where the pupil is non null
    protected double[] rho; // pupil modulus based on Zernike polynomials
    protected double[] phi; // pupil phase based on Zernike polynomials
    protected double[] psi; // defocus function
    protected Array4D cpxPsf; // fourier transform of the pupil function

    protected Shape cpxPsfShape;
    protected Shape aShape;
    protected Shape psf2DShape;

    protected  Object FFT2D;

    protected int nModulus;
    protected int nDefocus;
    protected int nPhase;

    private boolean para=true;

    public WideFieldModel(Shape psfShape,
            double NA, double lambda, double ni, double dxy, double dz, boolean radial, boolean single){
        this( psfShape,0, 0,
                NA,  lambda,  ni,  dxy,  dz,  radial,  single) ;
    }

    public WideFieldModel(Shape psfShape,int nPhase, int nModulus,
            double NA, double lambda, double ni, double dxy, double dz, boolean radial, boolean single) {
        // TODO Check size of input vectors.

        super(psfShape,  NA,  ni, dxy, dz,  radial, single);
        if(Nx != Ny){
            throw new IllegalArgumentException("Nx should equal Ny");
        }
        this.lambda = lambda;
        this.ni = ni;
        this.Nzern = 4;
        this.radius = NA/lambda;
        this.lambda_ni = ni/lambda;
        this.phi = new double[Ny*Nx];
        this.psi = new double[Ny*Nx];
        cpxPsfShape = new Shape(2,Nx, Ny,  Nz);
        aShape = new Shape(2,Nx, Ny);
        psf2DShape = new Shape(Nx, Ny);

        computeMaskPupil();
        computeZernike();

        if(single){
            FFT2D = new FloatFFT_2D(Nx, Ny);
        }else{
            FFT2D = new DoubleFFT_2D(Nx, Ny);
        }

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

    @SuppressWarnings("unused")
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
            ConcurrencyUtils.setNumberOfThreads(1);
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
                        ((FloatFFT_2D) FFT2D).complexForward(A);

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
                ConcurrencyUtils.setNumberOfThreads(1);
                for ( int iz = 0; iz < Nz; iz++)
                {
                    final int iz1 = iz;
                    Callable<GetPsfParaOut> callable = new Callable<GetPsfParaOut>() {
                        @Override
                        public GetPsfParaOut call() throws Exception {
                            GetPsfParaOut output = new GetPsfParaOut(Npix,iz1,single);

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
                            ((DoubleFFT_2D) FFT2D).complexForward(A);

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

                //    List<Output> outputs = new ArrayList<Output>();
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
                    ((DoubleFFT_2D) FFT2D).complexForward(A);


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
    public  DoubleShapedVector apply_J_modulus(DoubleShapedVector q)
    {
        int Ci;
        final int Npix = Nx*Ny;
        double defoc_scale = 0.;
        final double PSFNorm = 1.0/(Nx*Ny*Nz);
        double Aq[] = new double[2*Npix];
        double J[] = new double[Ny*Nx];
        final double NBeta =1./modulus_coefs.norm2();
        Array1D JRho = null;

        if(single){

        }else{
            JRho =  Double1D.create(modulusSpace.getShape());

            ((Double1D) JRho).fill(0.);
            if(para){
                //  final Double3D q1 = ((Double3D) q.asShapedArray());
                int threads = Runtime.getRuntime().availableProcessors();
                ExecutorService service = Executors.newFixedThreadPool(threads);

                List<Future<ApplyJPhaOut>> futures = new ArrayList<Future<ApplyJPhaOut>>();
                ConcurrencyUtils.setNumberOfThreads(1);
                for ( int iz = 0; iz < Nz; iz++)
                {

                    final Double2D qz = ((Double3D) q.asShapedArray()).slice(iz);
                    final int iz1 = iz;
                    Callable<ApplyJPhaOut> callable = new Callable<ApplyJPhaOut>() {
                        @Override
                        public ApplyJPhaOut call() throws Exception {


                            double defoc_scale=0;
                            double Aq[] = new double[2*Npix];
                            double defoc;
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
                            ((DoubleFFT_2D) FFT2D).complexForward(Aq);


                            for (int j = 0; j < Ny; j++)
                            {
                                for (int i = 0; i < Nx; i++)
                                {
                                    int in = i + j*Nx;
                                    if(maskPupil[in] == 1)
                                    {
                                        double idef= 1./psi[in];
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
                            ((Double1D) JRho).set(k,((Double1D) JRho).get(k)+ pout.grd[k]);
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
                DoubleFFT_2D FFT2D = new DoubleFFT_2D(Ny, Nx);

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
                    /*      for (int in = 0; in < Npix; in++)
                {
                    Ci = iz*Npix + in;
                    //   Aq[2*in] = cpxPsf[2*Ci]*q[Ci];
                    //   Aq[2*in + 1] = cpxPsf[2*Ci + 1]*q[Ci];
                    Double3D cpxPsfZ = (Double3D) cpxPsf.slice(iz);
                    Double1D qZ = ((Double3D) q.asShapedArray()).slice(iz).as1D();
                    Aq[2*in] =   cpxPsfZ.slice(0, 0).as1D().get(in)*qZ.get(in);
                    Aq[2*in + 1] = cpxPsfZ.slice(1, 0).as1D().get(in)*qZ.get(in);

                }*/

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
                    ((Double1D) JRho).set(k,2*PSFNorm*tmp*(1 - Math.pow(modulus_coefs.norm2()*NBeta,2))*NBeta);
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
    public DoubleShapedVector apply_J_phi(DoubleShapedVector q)
    {
        int Ci;
        final int Npix = Nx*Ny;
        final double PSFNorm = 1.0/(Nx*Ny*Nz);
        Double1D JPhi =  Double1D.create(phaseSpace.getShape());
        JPhi.fill(0.);
        if(para){
            //  final Double3D q1 = ((Double3D) q.asShapedArray());
            int threads = Runtime.getRuntime().availableProcessors();
            ExecutorService service = Executors.newFixedThreadPool(threads);

            List<Future<ApplyJPhaOut>> futures = new ArrayList<Future<ApplyJPhaOut>>();
            ConcurrencyUtils.setNumberOfThreads(1);
            for ( int iz = 0; iz < Nz; iz++)
            {

                final Double2D qz = ((Double3D) q.asShapedArray()).slice(iz);
                final int iz1 = iz;
                Callable<ApplyJPhaOut> callable = new Callable<ApplyJPhaOut>() {
                    @Override
                    public ApplyJPhaOut call() throws Exception {


                        double defoc_scale=0;
                        double Aq[] = new double[2*Npix];
                        double defoc;
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
                        ((DoubleFFT_2D) FFT2D).complexForward(Aq);


                        for (int j = 0; j < Ny; j++)
                        {
                            for (int i = 0; i < Nx; i++)
                            {
                                int in = i + j*Nx;
                                if(maskPupil[in] == 1)
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

            //    List<Output> outputs = new ArrayList<Output>();
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
                //            for (int in = 0; in < Npix; in++)
                //            {
                //                // Ci = iz*Npix + in;
                //                // Aq[2*in] = cpxPsf[2*Ci]*q[Ci];
                //                //Aq[2*in + 1] = cpxPsf[2*Ci + 1]*q[Ci];
                //                Double3D cpxPsfZ = (Double3D) cpxPsf.slice(iz);
                //                Double1D qZ = ((Double3D) q.asShapedArray()).slice(iz).as1D();
                //                Aq[2*in] =   cpxPsfZ.slice(0, 0).as1D().get(in)*qZ.get(in);
                //                Aq[2*in + 1] = cpxPsfZ.slice(1, 0).as1D().get(in)*qZ.get(in);
                //
                //            }

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
        return phaseSpace.create(JPhi );
    }


    /**
     * Apply the Jacobian matrix to go from  the PSF space to defocus coefficients space.
     * @param q : the gradient of some criterion in the PSF space
     * @return the gradient of this criterion in the defocus coefficients space.
     */
    @Override
    public DoubleShapedVector apply_J_defocus(DoubleShapedVector q)
    {

        // long startTime = System.currentTimeMillis();

        double scale_x = 1/(Nx*dxy);
        double scale_y = 1/(Ny*dxy);
        double defoc, tmpvar, idef, d0 = 0, d1 = 0, d2 = 0;
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
        if(para){
            //  final Double3D q1 = ((Double3D) q.asShapedArray());
            int threads = Runtime.getRuntime().availableProcessors();
            ExecutorService service = Executors.newFixedThreadPool(threads);

            List<Future<ApplyJDefOut>> futures = new ArrayList<Future<ApplyJDefOut>>();
            ConcurrencyUtils.setNumberOfThreads(1);
            for ( int iz = 0; iz < Nz; iz++)
            {

                final Double2D qz = ((Double3D) q.asShapedArray()).slice(iz);
                final int iz1 = iz;
                Callable<ApplyJDefOut> callable = new Callable<ApplyJDefOut>() {
                    @Override
                    public ApplyJDefOut call() throws Exception {


                        double defoc_scale=0;
                        double Aq[] = new double[2*Npix];
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
                                //double qin =  q1.get(ix, iy, iz1);
                                double qin =  qz.get(ix, iy);
                                Aq[2*in]=  ((Double4D) cpxPsf).get(0, ix, iy, iz1 )*qin;
                                Aq[2*in+1]=  ((Double4D) cpxPsf).get(1, ix, iy, iz1 )*qin;
                            }

                        }
                        /* Fourier transform of the pupil function A(z) */
                        ((DoubleFFT_2D) FFT2D).complexForward(Aq);
                        for (int j = 0; j < Ny; j++)
                        {
                            for (int i = 0; i < Nx; i++)
                            {
                                int in = i + j*Nx;
                                if(maskPupil[in] == 1)
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

            //    List<Output> outputs = new ArrayList<Output>();
            for (Future<ApplyJDefOut> future : futures) {
                ApplyJDefOut output;
                try {
                    output = future.get();
                    //  System.out.println("d0: "+d0 + " d1: " + d1 + " d2: "+d2);
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


            double Aq[] = new double[2*Npix];
            DoubleFFT_2D FFT2D = new DoubleFFT_2D(Ny, Nx);

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
                        double qin =  ((Double3D) q.asShapedArray()).get(ix, iy, iz);
                        Aq[2*in]=  ((Double4D) cpxPsf).get(0, ix, iy, iz )*qin;
                        Aq[2*in+1]=  ((Double4D) cpxPsf).get(1, ix, iy, iz )*qin;
                    }

                }


                FFT2D.complexForward(Aq);


                for (int j = 0; j < Ny; j++)
                {
                    for (int i = 0; i < Nx; i++)
                    {
                        int in = i + j*Nx;
                        if(maskPupil[in] == 1)
                        {
                            // Ci = iz*Npix + in;
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


        //        long stopTime = System.currentTimeMillis();
        //        long elapsedTime = stopTime - startTime;
        //        System.out.println("time: "+elapsedTime);

        return  defocusSpace.create(Double1D.wrap(grd, defocusSpace.getShape()));

    }


    /** Determine the map where the pupil in non null.  It sets  maskPupil and its area pupil_area
     */
    private void computeMaskPupil()
    {
        maskPupil = new double[Nx*Ny];
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
                    maskPupil[nx + ny*Nx] = 1;
                    pupil_area += 1;

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
            Arrays.fill(grd, 0);
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
                if (maskPupil[nxy] == 1)
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
                        maskPupil[nxy] = 0;
                    }
                    else
                    {
                        psi[nxy] = Math.sqrt(q);
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
            // TODO Error
        }

        int Npix = Nx*Ny;
        rho = new double[Npix];
        //    double betaNorm = 1./(Math.sqrt(MathUtils.innerProd(modulus_coefs, modulus_coefs)));
        double betaNorm = 1./( beta.norm2());
        for(int in = 0; in < Npix; in++)
        {
            if (maskPupil[in] == 1)
            {

                for (int n = 0; n < beta.getNumber(); ++n)
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
            if (maskPupil[in] == 1)
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
    public double[] getMaskPupil() {
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

        System.out.println("----PUPIL----");
        MathUtils.stat(maskPupil);
        System.out.println();

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



}

