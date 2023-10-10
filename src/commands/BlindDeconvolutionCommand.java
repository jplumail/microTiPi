package commands;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;
import java.util.Locale;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import microTiPi.epifluorescence.WideFieldModel;
import microTiPi.microUtils.BlindDeconvJob;
import microTiPi.microscopy.PSF_Estimation;
import mitiv.array.ShapedArray;
import mitiv.conv.WeightedConvolutionCost;
import mitiv.cost.DifferentiableCostFunction;
import mitiv.cost.HyperbolicTotalVariation;
import mitiv.io.DataFormat;
import mitiv.jobs.DeconvolutionJob;
import mitiv.linalg.shaped.DoubleShapedVectorSpace;
import mitiv.linalg.shaped.FloatShapedVectorSpace;
import mitiv.linalg.shaped.ShapedVectorSpace;
import mitiv.weights.weightsFromModel;

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

    @Option(name = "-lambda", usage = "Wavelength in m")
    private double lambda = 0.00000054; // 540nm

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


    public static void main(String[] args) {

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
        

        ShapedArray dataArray = DataFormat.load(inputName);
        job.stream.format("dataArray shape: %s\n", dataArray.getShape());
        WideFieldModel pupil = new WideFieldModel(dataArray.getShape(), job.nPhase, job.nModulus, job.NA, job.lambda, job.ni, job.dxy, job.dz, job.radial, job.single);
        PSF_Estimation psfEstimation = new PSF_Estimation(pupil);

        ShapedVectorSpace dataSpace, objectSpace;
        if (job.single) {
            objectSpace = new FloatShapedVectorSpace(dataArray.getShape());
            dataSpace = new FloatShapedVectorSpace(dataArray.getShape());
        }
        else {
            objectSpace = new DoubleShapedVectorSpace(dataArray.getShape());
            dataSpace = new DoubleShapedVectorSpace(dataArray.getShape());
        }
        double[] scale = {1, 1, job.dz / job.dxy};
        DifferentiableCostFunction fprior = new HyperbolicTotalVariation(objectSpace, job.epsilon, scale);
        WeightedConvolutionCost fdata =  WeightedConvolutionCost.build(objectSpace, dataSpace);
        DeconvolutionJob deconvolver = new DeconvolutionJob(fdata, job.mu, fprior, !job.negativity, job.nbIterDeconv, null, null); // hooks to null

        weightsFromModel wghtUpdt = new weightsFromModel(dataArray, null); // badPixArray to null, no pb?

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
