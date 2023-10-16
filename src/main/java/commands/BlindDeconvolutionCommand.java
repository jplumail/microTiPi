package commands;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.ByteBuffer;
import java.util.List;
import java.util.Locale;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.mitiv.microTiPi.epifluorescence.WideFieldModel;
import org.mitiv.microTiPi.microUtils.BlindDeconvJob;
import org.mitiv.microTiPi.microscopy.PSF_Estimation;
import org.mitiv.TiPi.array.Array3D;
import org.mitiv.TiPi.array.ArrayFactory;
import org.mitiv.TiPi.array.ArrayUtils;
import org.mitiv.TiPi.array.DoubleArray;
import org.mitiv.TiPi.array.FloatArray;
import org.mitiv.TiPi.array.ShapedArray;
import org.mitiv.TiPi.base.Shape;
import org.mitiv.TiPi.conv.WeightedConvolutionCost;
import org.mitiv.TiPi.cost.DifferentiableCostFunction;
import org.mitiv.TiPi.cost.HyperbolicTotalVariation;
import org.mitiv.TiPi.jobs.DeconvolutionJob;
import org.mitiv.TiPi.linalg.shaped.DoubleShapedVectorSpace;
import org.mitiv.TiPi.linalg.shaped.FloatShapedVectorSpace;
import org.mitiv.TiPi.linalg.shaped.ShapedVectorSpace;
import org.mitiv.TiPi.utils.FFTUtils;
import org.mitiv.io.DeconvHook;
import org.mitiv.io.NullImager;
import org.mitiv.TiPi.weights.WeightFactory;

import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.ImageWriter;
import loci.formats.in.OMETiffReader;
import loci.formats.meta.IMetadata;
import loci.formats.out.OMETiffWriter;
import loci.formats.services.OMEXMLService;
import ome.xml.model.enums.DimensionOrder;
import ome.xml.model.enums.PixelType;
import ome.xml.model.primitives.PositiveInteger;

public class BlindDeconvolutionCommand {

    @Option(name = "-nbloops", usage = "number of loops")
    private int loops = 1;

    @Option(name = "-maxIterDefocus", usage = "Max number of iterations for defocus.")
    private int maxIterDefocus = 20;

    @Option(name = "-maxIterPhase", usage = "Max number of iterations for phase.")
    private int maxIterPhase = 20;

    @Option(name = "-maxIterModulus", usage = "Max number of iterations for modulus.")
    private int maxIterModulus = 0;

    // WidefieldModel args
    @Option(name = "-nPhase", usage = "Number of zernike describing the pupil phase.")
    private int nPhase = 19;

    @Option(name = "-nModulus", usage = "Number of zernike describing the pupil modulus.")
    private int nModulus = 0;

    @Option(name = "-NA", usage = "Numerical aperture")
    private double NA = 1.4;

    @Option(name = "-lambda", usage = "Wavelength in nm")
    private double lambda = 540; // 540nm

    @Option(name = "-ni", usage = "Refractive index")
    private double ni = 1.518;

    @Option(name = "-dxy", usage = "Lateral pixel size in nm")
    private double dxy = 1.;

    @Option(name = "-dz", usage = "Axial pixel size in nm")
    private double dz = 1.;

    @Option(name = "-radial", usage = "Radial option")
    private boolean radial;

    @Option(name = "-single", usage = "Compute in single precision")
    private boolean single;

    // Deconvolution args
    @Option(name = "-epsilon", usage = "threshold of Hyperbolic TV")
    private double epsilon = 0.01;

    @Option(name = "-mu", usage = "Value of mu. Must be positive")
    private double mu = 1.;

    @Option(name = "-negativity", usage = "Allow negativity")
    private boolean negativity;

    @Option(name = "-nbIterDeconv", usage = "Number of iterations for deconvolution")
    private int nbIterDeconv = 50;

    // Misc
    @Option(name = "-debug", usage = "debug flag")
    private boolean debug;

    @Option(name = "-help", aliases = {"--help", "-h", "-?"}, usage = "Display help.")
    private boolean help;

    @Argument
    private List<String> arguments;

    static private void usage(CmdLineParser parser, int code) {
        PrintStream stream = (code == 0 ? System.out : System.err);
        stream.println("Usage: blinddeconv [OPTIONS] INPUT OUTPUT PSF");
        if (code == 0) {
            stream.println("Options:");
            parser.getProperties().withUsageWidth(80);
            parser.printUsage(stream);
        } else {
            stream.println("Try option -help for a more complete description of options.");
        }
        System.exit(code);
    }

    private static WideFieldModel buildPupil(BlindDeconvolutionCommand job, Shape psfShape) {
        WideFieldModel pupil = new WideFieldModel(psfShape, job.nPhase, job.nModulus, job.NA, job.lambda*1E-9, job.ni, job.dxy*1E-9, job.dz*1E-9, job.radial, job.single);
        return pupil;
    }


    public static void main(String[] args) throws FormatException, IOException, DependencyException, ServiceException {

        // Switch to "US" locale to avoid problems with number formats.
        Locale.setDefault(Locale.US);

        // Parse options.
        BlindDeconvolutionCommand job = new BlindDeconvolutionCommand();
        CmdLineParser parser = new CmdLineParser(job);
        try {
            parser.parseArgument(args);
        } catch (CmdLineException e) {
            System.err.format("Error: %s\n", e.getMessage());
            usage(parser, 1);
        }
        if (job.help) {
            usage(parser, 0);
        }

        // Deal with remaining arguments.
        int size = (job.arguments == null ? 0 : job.arguments.size());
        if (size != 3) {
            System.err.format("Too %s arguments.\n", (size < 2 ? "few" : "many"));
            usage(parser, 1);
        }
        String inputName = job.arguments.get(0);
        String outputName = job.arguments.get(1);
        String psfName = job.arguments.get(2);
        
        // Read input file
        ShapedArray dataArray = readOMETiff(inputName, job.single);

        // Compute output shape
        ShapedVectorSpace dataSpace, objectSpace;
        int Nx, Ny, Nz;
        Nx = FFTUtils.bestDimension(dataArray.getDimension(0));
        Ny = FFTUtils.bestDimension(dataArray.getDimension(1));
        Nz = FFTUtils.bestDimension(dataArray.getDimension(2));
        Shape outputShape = new Shape(Nx, Ny, Nz);

        // build objArray, the deconvolved image
        ShapedArray objArray = dataArray.copy();
        objArray = ArrayUtils.extract(objArray, outputShape, 0.); //Padding to the right size

        // build pupil, psf
        WideFieldModel pupil = buildPupil(job, outputShape);
        PSF_Estimation psfEstimation = new PSF_Estimation(pupil);
        ShapedArray psfArray = ArrayUtils.roll( pupil.getPsf() );

        // build hooks
        NullImager imager = new NullImager();
        DeconvHook dHook = new DeconvHook(imager, outputShape,null, job.debug);
        DeconvHook dHookfinal = new DeconvHook(imager, outputShape,"Deconvolved", job.debug);

        // build deconvolver
        if (job.single) {
            dataSpace = new FloatShapedVectorSpace(dataArray.getShape());
            objectSpace = new FloatShapedVectorSpace(outputShape);
        }
        else {
            dataSpace = new DoubleShapedVectorSpace(dataArray.getShape());
            objectSpace = new DoubleShapedVectorSpace(outputShape);
        }
        double[] scale = {1, 1, job.dz / job.dxy};
        DifferentiableCostFunction fprior = new HyperbolicTotalVariation(objectSpace, job.epsilon, scale);
        WeightedConvolutionCost fdata =  WeightedConvolutionCost.build(objectSpace, dataSpace);
        fdata.setData(dataArray);
        fdata.setPSF(psfArray);
        DeconvolutionJob deconvolver = new DeconvolutionJob(fdata, job.mu, fprior, !job.negativity, job.nbIterDeconv, dHook, dHookfinal); // hooks to null

        // build wgtArray, weights
        // Compute variance method, TODO: implement others
        double gamma = 1.;
        double sigma = 10.;
        double alpha = gamma;
        double beta = (sigma/gamma)*(sigma/gamma);
        ShapedArray wgtArray = WeightFactory.computeWeightsFromData(dataArray, alpha, beta);
        WeightFactory.normalize(wgtArray);

        // update weights
        fdata.setWeights(wgtArray,true);
        psfEstimation.setWeight(ArrayUtils.pad(wgtArray,outputShape));

        // update psfEstimation
        psfEstimation.setData(objArray);
        psfEstimation.enablePositivity(false);
        psfEstimation.setAbsoluteTolerance(0.0);

        // build bdec
        BlindDeconvJob bdec = new BlindDeconvJob(
            job.loops,
            pupil.getParametersFlags(),
            new int[] {job.maxIterDefocus, job.maxIterPhase, job.maxIterModulus},
            psfEstimation,
            deconvolver,
            null, // TODO: wghtUpdt to null because Compute variance method
            job.debug
        );

        // run blinddeconv
        // System.out.println("Blind deconv");
        // long startTime = System.nanoTime();
        objArray = bdec.blindDeconv(objArray);
        // long endTime = System.nanoTime();
        // long elapsedTime = endTime - startTime;
        // double elapsedTimeInSeconds = (double) elapsedTime / 1_000_000_000.0;
        // System.out.println("Elapsed Time: " + elapsedTimeInSeconds + " seconds");
        
        // get the psf
        pupil = ((WideFieldModel) psfEstimation.getModel());
        pupil.freeMem();

        // save arrays
        saveArray(outputName, objArray, job.single);
        saveArray(psfName, ArrayUtils.roll(pupil.getPsf()), job.single);
    }

    private static ShapedArray readOMETiff(String path, boolean single) throws FormatException, IOException {
        OMETiffReader reader = new OMETiffReader();
        reader.setId(path);
        if (reader.getSeriesCount()>1 || reader.getSizeT()>1 || reader.getSizeC()>1) {
            reader.close();
            throw new FormatException("File no good shape (Series:%d, T:%d, C:%d)".formatted(reader.getSeriesCount(), reader.getSizeT(), reader.getSizeC()));
        }
        reader.setSeries(0);
        int bitsPerPixel = reader.getBitsPerPixel();
        int sizeX = reader.getSizeX();
        int sizeY = reader.getSizeY();
        int sizeZ = reader.getSizeZ();
        // Calculate the size in bits
        int bufferSizeInBits = bitsPerPixel * sizeX * sizeY * sizeZ;

        // Allocate the ByteBuffer
        ByteBuffer buffer = ByteBuffer.allocate(bufferSizeInBits);
        for (int i=0; i<reader.getSizeZ(); i++) {
            byte[] plane = reader.openBytes(i);
            buffer.put(plane);
        }
        ShapedArray dataArray;
        if (single) {
            ShapedArray shapedArray = ArrayFactory.wrap(buffer.array(), reader.getSizeX(), reader.getSizeY(), reader.getSizeZ());
            FloatArray floatArray = shapedArray.toFloat();
            floatArray.scale(1/floatArray.max());
            dataArray = (ShapedArray) floatArray;
        } else {
            ShapedArray shapedArray = ArrayFactory.wrap(buffer.array(), reader.getSizeX(), reader.getSizeY(), reader.getSizeZ());
            DoubleArray doubleArray = shapedArray.toDouble();
            doubleArray.scale(1/doubleArray.max());
            dataArray = (ShapedArray) doubleArray;
        }
        reader.close();
        return dataArray;
    }

    private static void saveArray(String path, ShapedArray arr, boolean single)
    throws DependencyException, ServiceException, FormatException, IOException {
        // System.out.println("Start saving");
        // long startTime = System.nanoTime();
        ServiceFactory factory = new ServiceFactory();
        OMEXMLService service = factory.getInstance(OMEXMLService.class);
        IMetadata omexml = service.createOMEXMLMetadata();
        omexml.setImageID("Image:0", 0);
        omexml.setPixelsID("Pixels:0", 0);
        omexml.setPixelsBinDataBigEndian(Boolean.FALSE, 0, 0);
        omexml.setPixelsDimensionOrder(DimensionOrder.XYCZT, 0);
        if (single) {
            omexml.setPixelsType(PixelType.FLOAT, 0);
        } else {
            omexml.setPixelsType(PixelType.DOUBLE, 0);
        }
        omexml.setPixelsSizeX(new PositiveInteger(arr.getDimension(0)), 0);
        omexml.setPixelsSizeY(new PositiveInteger(arr.getDimension(1)), 0);
        omexml.setPixelsSizeZ(new PositiveInteger(arr.getDimension(2)), 0);
        omexml.setPixelsSizeT(new PositiveInteger(1), 0);
        omexml.setPixelsSizeC(new PositiveInteger(1), 0);
        omexml.setChannelID("Channel:0:0", 0, 0);
        omexml.setChannelSamplesPerPixel(new PositiveInteger(1),0, 0);

        ImageWriter imwriter = new ImageWriter();
        imwriter.setMetadataRetrieve(omexml);
        imwriter.setId(path);
        OMETiffWriter writer = (OMETiffWriter) imwriter.getWriter();
        // tiffWriter.setCompression(OMETiffWriter.COMPRESSION_UNCOMPRESSED);
        Array3D data = (Array3D) arr;
        // long rowsPerStrip = 64000; // use all rows
        // long[] rowsPerStripArray = new long[]{ rowsPerStrip };
        if (single){
            for (int image=0; image<arr.getDimension(2); image++) {
                // IFD ifd = new IFD();
                //ifd.put( IFD.ROWS_PER_STRIP, rowsPerStripArray );
                float[] plane = (float[]) data.slice(image, 2).flatten(true);
                ByteBuffer bb = ByteBuffer.allocate(plane.length * 4);
                for (float d: plane) {
                    bb.putFloat(d);
                }
                writer.saveBytes(image, bb.array());
            }
        } else{
            for (int image=0; image<arr.getDimension(2); image++) {
                // IFD ifd = new IFD();
                //ifd.put( IFD.ROWS_PER_STRIP, rowsPerStripArray );
                double[] plane = (double[]) data.slice(image, 2).flatten(true);
                ByteBuffer bb = ByteBuffer.allocate(plane.length * 8);
                for (double d: plane) {
                    bb.putDouble(d);
                }
                writer.saveBytes(image, bb.array());
            }
        }
        writer.close();
        imwriter.close();
        // long endTime = System.nanoTime();
        // long elapsedTime = endTime - startTime;
        // double elapsedTimeInSeconds = (double) elapsedTime / 1_000_000_000.0;
        // System.out.println("Elapsed Time: " + elapsedTimeInSeconds + " seconds");

    }
}
