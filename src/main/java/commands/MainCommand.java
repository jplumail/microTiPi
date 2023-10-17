package commands;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.ByteBuffer;
import java.util.Arrays;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.mitiv.TiPi.array.Array3D;
import org.mitiv.TiPi.array.ArrayFactory;
import org.mitiv.TiPi.array.DoubleArray;
import org.mitiv.TiPi.array.FloatArray;
import org.mitiv.TiPi.array.ShapedArray;
import org.mitiv.TiPi.io.ColorModel;

import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.IFormatWriter;
import loci.formats.ImageReader;
import loci.formats.ImageWriter;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import ome.xml.model.enums.DimensionOrder;
import ome.xml.model.enums.PixelType;
import ome.xml.model.primitives.PositiveInteger;

public class MainCommand {
    private PrintStream stream = System.out;

    @Option(name="prog", usage="choose the program to use: deconv or blinddeconv")
    private String arg1;

    @Argument
    private String args;

    public static void main(String[] args) throws FormatException, IOException, DependencyException, ServiceException {
        MainCommand job = new MainCommand();
        if (args.length > 1){
            String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
            if (args[0].equals("deconv")) {
                EdgePreservingDeconvolutionCommand.main(newArgs);
                return;
            } else if (args[0].equals("blinddeconv")) {
                BlindDeconvolutionCommand.main(newArgs);
                return;
            }
        }
        CmdLineParser parser = new CmdLineParser(job);
        job.stream.println("Usage: microtipi prog [OPTIONS] INPUT OUTPUT");
        parser.printUsage(job.stream);
    }

    static ShapedArray readOMETiffToArray(String path, boolean single) throws FormatException, IOException {
        ImageReader reader = new ImageReader();
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

    static void saveArrayToOMETiff(String path, ShapedArray arr, boolean single)
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
        IFormatWriter writer = imwriter.getWriter();
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

    public static ShapedArray loadData(String name, boolean single) throws FormatException, IOException {
        ShapedArray arr = readOMETiffToArray(name, single);
        ColorModel colorModel = ColorModel.guessColorModel(arr);
        if (colorModel == ColorModel.NONE) {
            return (single ? arr.toFloat() :  arr.toDouble());
        } else {
            return (single
                    ? ColorModel.filterImageAsFloat(arr, ColorModel.GRAY)
                            : ColorModel.filterImageAsDouble(arr, ColorModel.GRAY));
        }
    }
}
