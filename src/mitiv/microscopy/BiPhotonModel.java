package mitiv.microscopy;

public class BiPhotonModel extends MicroscopeModel {
    Widefi

    public BiPhotonModel(double NA, double lambda, double ni, double dxy, double dz, int Nx, int Ny, int Nz,
            boolean radial) {
        super(NA, lambda, ni, dxy, dz, Nx, Ny, Nz, radial);
    }

    @Override
    void computePSF() {
        // TODO Auto-generated method stub

    }
}

