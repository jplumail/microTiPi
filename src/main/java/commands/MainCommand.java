package commands;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;

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
}
