package commands;

import java.io.FileNotFoundException;
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
import org.mitiv.TiPi.array.ArrayFactory;
import org.mitiv.TiPi.array.ArrayUtils;
import org.mitiv.TiPi.array.ShapedArray;
import org.mitiv.TiPi.base.Shape;
import org.mitiv.TiPi.conv.WeightedConvolutionCost;
import org.mitiv.TiPi.cost.DifferentiableCostFunction;
import org.mitiv.TiPi.cost.HyperbolicTotalVariation;
import org.mitiv.TiPi.io.DataFormat;
import org.mitiv.TiPi.jobs.DeconvolutionJob;
import org.mitiv.TiPi.linalg.shaped.DoubleShapedVectorSpace;
import org.mitiv.TiPi.linalg.shaped.FloatShapedVectorSpace;
import org.mitiv.TiPi.linalg.shaped.ShapedVectorSpace;
import org.mitiv.TiPi.utils.FFTUtils;
import org.mitiv.TiPi.utils.HistoMap;
import org.mitiv.TiPi.weights.weightsFromModel;

import loci.formats.FormatException;
import loci.formats.in.OMETiffReader;

public class BlindDeconvolutionCommand {
    private PrintStream stream = System.out;

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
        stream.println("Usage: blinddeconv [OPTIONS] INPUT OUTPUT");
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


    public static void main(String[] args) throws FormatException, IOException {

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
        if (size != 2) {
            System.err.format("Too %s arguments.\n", (size < 2 ? "few" : "many"));
            usage(parser, 1);
        }
        String inputName = job.arguments.get(0);
        String outputName = job.arguments.get(1);
        
        try (OMETiffReader reader = new OMETiffReader()) {
            reader.setId(inputName);
            if (reader.getSeriesCount()>1 || reader.getSizeT()>1 || reader.getSizeC()>1) {
                throw new FormatException("File no good shape (Series:%d, T:%d, C:%d)".formatted(reader.getSeriesCount(), reader.getSizeT(), reader.getSizeC()));
            }
            reader.setSeries(0);
            int bitsPerPixel = reader.getBitsPerPixel(); // Replace 'reader' with your actual reader instance
            int sizeX = reader.getSizeX(); // Replace 'reader' with your actual reader instance
            int sizeY = reader.getSizeY(); // Replace 'reader' with your actual reader instance
            int sizeZ = reader.getSizeZ(); // Replace 'reader' with your actual reader instance
            // Calculate the size in bits
            int bufferSizeInBits = bitsPerPixel * sizeX * sizeY * sizeZ;

            // Allocate the ByteBuffer
            ByteBuffer buffer = ByteBuffer.allocate((int) (bufferSizeInBits)); // Convert bits to bytes
            for (int i=0; i<reader.getSizeZ(); i++) {
                byte[] plane = reader.openBytes(i);
                buffer.put(plane);
            }
            ShapedArray dataArray = ArrayFactory.wrap(buffer.array(), reader.getSizeX(), reader.getSizeY(), reader.getSizeZ());

            job.stream.format("dataArray shape: %s\n", dataArray.getShape());
            Shape psfShape = new Shape(64, 64, 64); // possible?
            WideFieldModel pupil = buildPupil(job, psfShape);
            
            PSF_Estimation psfEstimation = new PSF_Estimation(pupil);
            ShapedArray psfArray = ArrayUtils.roll( pupil.getPsf() );


            ShapedVectorSpace dataSpace, objectSpace;
            int Nx, Ny, Nz;
            Nx = FFTUtils.bestDimension(sizeX);
            Ny = FFTUtils.bestDimension(sizeY);
            Nz= FFTUtils.bestDimension(sizeZ);
            Shape outputShape = new Shape(Nx, Ny, Nz);
            if (job.single) {
                dataSpace = new FloatShapedVectorSpace(dataArray.getShape());
                objectSpace = new FloatShapedVectorSpace(outputShape);
            }
            else {
                dataSpace = new DoubleShapedVectorSpace(dataArray.getShape());
                objectSpace = new DoubleShapedVectorSpace(outputShape);
            }
            dataArray = ArrayUtils.extract(dataArray, outputShape, 0.);
            double[] scale = {1, 1, job.dz / job.dxy};
            DifferentiableCostFunction fprior = new HyperbolicTotalVariation(objectSpace, job.epsilon, scale);
            WeightedConvolutionCost fdata =  WeightedConvolutionCost.build(objectSpace, dataSpace);
            DeconvolutionJob deconvolver = new DeconvolutionJob(fdata, job.mu, fprior, !job.negativity, job.nbIterDeconv, null, null); // hooks to null

            weightsFromModel wghtUpdt = new weightsFromModel(dataArray, null); // badPixArray to null
            HistoMap hm = new HistoMap(modelArray, dataArray, null);
            gain.setValue(hm.getAlpha());
            noise.setValue(Math.sqrt(hm.getBeta())/gain.getValue());
            wgtArray = hm.computeWeightMap(modelArray);
            wghtUpdt.update(deconvolver); // compute weights: option 5 Automatic variance estimation TODO: implement other methods (https://github.com/FerreolS/tipi4icy/blob/master/src/plugins/ferreol/demics/DEMICSPlug.java#L123)
            fdata.setData(dataArray);
            fdata.setPSF(psfArray);

            int[] maxiter = {job.maxIterDefocus, job.maxIterPhase, job.maxIterModulus};
            BlindDeconvJob bdec = new BlindDeconvJob(job.loops, pupil.getParametersFlags(), maxiter, psfEstimation, deconvolver, wghtUpdt, job.debug);

            dataArray = bdec.blindDeconv(dataArray);
            
            try {
                DataFormat.save(dataArray, outputName);
            } catch (FileNotFoundException e) {
                job.stream.format("File %s not found.\n", outputName);
            } catch (IOException e){
                job.stream.format("Could not write to file %s.\n", outputName);
            }
        }
    }
}
