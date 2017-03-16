package microTiPi.microscopy;

/**
 * Define a class for meta-data
 * @author ferreol
 *
 */
public class MicroscopeMetadata {
    //Just a public object with all psf values inside
    /**
     * Lateral size of a pixel in nm
     */
    public double dxy     = 64.5;
    /**
     * axial sampling step size
     */
    public double dz      = 160;
    /**
     * number of pixels along x or y
     */
    public int    nxy     = 256;
    /**
     * number of planes
     */
    public int    nz      = 128;
    /**
     * Numerical aperture
     */
    public double na      = 1.4;
    /**
     * Wavelength in nm
     */
    public double lambda  = 542;
    /**
     * Refractive index of the immersion medium
     */
    public double ni      = 1.518;

    @Override
    public String toString(){
        return new String("dxy: "+dxy+" dz: "+dz+" nxy: "+nxy+" nz: "+nz+" na "+na+" lambda "+lambda+" ni "+ni);
    }
}
