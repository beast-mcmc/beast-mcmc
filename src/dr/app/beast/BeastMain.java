/*
 * BeastMain.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.app.beast;

import beagle.BeagleFlag;
import beagle.BeagleInfo;
import dr.app.beauti.BeautiMenuBarFactory;
import dr.app.checkpoint.BeastCheckpointer;
import dr.app.plugin.Plugin;
import dr.app.plugin.PluginLoader;
import dr.app.util.Arguments;
import dr.app.util.Utils;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmcmc.MCMCMC;
import dr.inference.mcmcmc.MCMCMCOptions;
import dr.inference.operators.OperatorSchedule;
import dr.math.MathUtils;
import dr.util.*;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParser;
import jam.util.IconUtils;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.List;
import java.util.logging.*;

public class BeastMain {

    private final static Version version = new BeastVersion();

    public static final double DEFAULT_DELTA = 1.0;
    public static final int DEFAULT_SWAP_CHAIN_EVERY = 100;

    private static final String CITATION_FILE_SUFFIX = ".citations.txt";

    static class BeastConsoleApp extends jam.console.ConsoleApplication {
        XMLParser parser = null;
        public final static String backgroundColor =  "#35484F";
        public final static String foregroundColor = "#CBB944";

        public BeastConsoleApp(String nameString, String titleString, String aboutString, javax.swing.Icon icon) throws IOException {
            super(nameString, titleString, aboutString, Color.decode(backgroundColor), Color.decode(foregroundColor), icon, false);
            getDefaultFrame().setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
        }

        public void doStop() {
            Iterator iter = parser.getThreads();
            while (iter.hasNext()) {
                Thread thread = (Thread) iter.next();
                //noinspection removal
                thread.stop(); // http://java.sun.com/j2se/1.5.0/docs/guide/misc/threadPrimitiveDeprecation.html
            }
        }

        public void setTitle(String title) {
            getDefaultFrame().setTitle(title);
        }
    }

    public BeastMain(File inputFile, BeastConsoleApp consoleApp, int maxErrorCount, final boolean verbose,
                     boolean parserWarning, boolean strictXML, List<String> additionalParsers,
                     MCMCMCOptions mc3Options) {


        if (inputFile == null) {
            throw new RuntimeException("Error: no input file specified");
        }

        String fileName = inputFile.getName();

        final Logger infoLogger = Logger.getLogger("dr.apps.beast");
        try {

            FileReader fileReader = new FileReader(inputFile);

            XMLParser parser = new BeastParser(new String[]{fileName}, additionalParsers, verbose, parserWarning, strictXML, version);

            if (consoleApp != null) {
                consoleApp.parser = parser;
            }

            // Add a handler to handle warnings and errors. This is a ConsoleHandler
            // so the messages will go to StdOut..
            Logger logger = Logger.getLogger("dr");

            Handler messageHandler = new MessageLogHandler();
            messageHandler.setFilter(new Filter() {
                public boolean isLoggable(LogRecord record) {
                    return record.getLevel().intValue() < Level.WARNING.intValue();
                }
            });
            logger.addHandler(messageHandler);

            // Add a handler to handle warnings and errors. This is a ConsoleHandler
            // so the messages will go to StdErr..
            Handler errorHandler = new ConsoleHandler();
            errorHandler.setFilter(new Filter() {
                public boolean isLoggable(LogRecord record) {
                    if (verbose) {
                        return record.getLevel().intValue() >= Level.WARNING.intValue();
                    } else {
                        return record.getLevel().intValue() >= Level.SEVERE.intValue();
                    }
                }
            });
            infoLogger.addHandler(errorHandler);

            if (Boolean.parseBoolean(System.getProperty("output_citations"))) {
                String citationFileName;
                if (System.getProperty("citations.filename") != null) {
                    citationFileName = System.getProperty("citations.filename");
                } else {
                    citationFileName = fileName.substring(0, fileName.toLowerCase().indexOf(".xml")) + CITATION_FILE_SUFFIX;
                }
                FileOutputStream citationStream = new FileOutputStream(FileHelpers.getFile(citationFileName));
                //Handler citationHandler = new MessageLogHandler(citationStream);
                Handler citationHandler = CitationLogHandler.getHandler(citationStream);
                //Logger.getLogger("dr.app.beast").addHandler(citationHandler);
                Logger.getLogger("dr.util").addHandler(citationHandler);
            }

            logger.setUseParentHandlers(false);

            infoLogger.info("Parsing XML file: " + fileName);
            infoLogger.info("  File encoding: " + fileReader.getEncoding());

            // This is a special logger that is for logging numerical and statistical errors
            // during the MCMC run. It will tolerate up to maxErrorCount before throwing a
            // RuntimeException to shut down the run.
            Logger errorLogger = Logger.getLogger("error");
            messageHandler = new ErrorLogHandler(maxErrorCount);
            messageHandler.setLevel(Level.WARNING);
            errorLogger.addHandler(messageHandler);

            PluginLoader.loadPlugins(parser);

            // Install the checkpointer. This creates a factory that returns
            // appropriate savers and loaders according to the user's options.
            //new BeastCheckpointer();

            // if any of these properties is not null, overrule XML checkpointing settings
            if (System.getProperty("save.state.file", null) != null ||
                    System.getProperty("save.state.at", null) != null ||
                    System.getProperty("save.state.every", null) != null ||
                    System.getProperty("save.state.time", null) != null) {
                System.setProperty("checkpointOverrule","true");
            }

            if (Boolean.parseBoolean(System.getProperty("checkpointOverrule", "false"))) {
                BeastCheckpointer.getInstance(null, -1, -1, false);
                Logger.getLogger("dr.apps.beast").info("Overriding checkpointing settings in the provided XML file");
            }

            if (mc3Options == null) {

                // just parse the file running all threads...
                parser.parse(fileReader, true);

            } else {

                double[] chainTemperatures = mc3Options.getChainTemperatures();

                int chainCount = chainTemperatures.length;
                MCMC[] chains = new MCMC[chainCount];
//                MCMCMCOptions options = new MCMCMCOptions(chainTemperatures, swapChainsEvery);

                Logger.getLogger("dr.apps.beast").info("Starting cold chain plus hot chains with temperatures: ");
                for (int i = 1; i < chainTemperatures.length; i++) {
                    Logger.getLogger("dr.apps.beast").info("Hot Chain " + i + ": " + chainTemperatures[i]);
                }

                Logger.getLogger("dr.apps.beast").info("Parsing XML file: " + fileName);

                // parse the file for the initial cold chain returning the MCMC object
                chains[0] = (MCMC) parser.parse(fileReader, MCMC.class);
                if (chains[0] == null) {
                    throw new dr.xml.XMLParseException("BEAST XML file is missing an MCMC element");
                }
                fileReader.close();

                chainTemperatures[0] = 1.0; // TODO Should perform in Mc3Options constructor

                for (int i = 1; i < chainCount; i++) {
                    // parse the file once for each hot chain
                    fileReader = new FileReader(inputFile);

                    // turn off all messages for subsequent reads of the file (they will be the same as the
                    // first time).
                    logger.setLevel(Level.OFF);

                    parser = new BeastParser(new String[]{fileName}, additionalParsers, verbose, parserWarning, strictXML, version);

                    PluginLoader.loadPlugins(parser);
                    
                    chains[i] = (MCMC) parser.parse(fileReader, MCMC.class);
                    if (chains[i] == null) {
                        throw new dr.xml.XMLParseException("BEAST XML file is missing an MCMC element");
                    }
                    fileReader.close();
                }

                // restart messages
                logger.setLevel(Level.ALL);

                MCMCMC mc3 = new MCMCMC(chains, mc3Options);
//                Thread thread = new Thread(mc3);
//                thread.start();
                mc3.run();
            }

        } catch (java.io.IOException ioe) {
            infoLogger.severe("File error: " + ioe.getMessage());
            throw new RuntimeException("Terminate");
        } catch (org.xml.sax.SAXParseException spe) {
            if (spe.getMessage() != null && spe.getMessage().equals("Content is not allowed in prolog")) {
                infoLogger.severe("Parsing error - the input file, " + fileName + ", is not a valid XML file.");
            } else {
                infoLogger.severe("Error running file: " + fileName);
                infoLogger.severe("Parsing error - poorly formed XML (possibly not an XML file):\n" +
                        spe.getMessage());
            }
            throw new RuntimeException("Terminate");
        } catch (org.w3c.dom.DOMException dome) {
            infoLogger.severe("Error running file: " + fileName);
            infoLogger.severe("Parsing error - poorly formed XML:\n" +
                    dome.getMessage());
            throw new RuntimeException("Terminate");
        } catch (dr.xml.XMLParseException pxe) {
            // Leave the printing of the stack trace until the end - too noisy otherwise
            //pxe.printStackTrace(System.err);
            if (pxe.getMessage() != null && pxe.getMessage().equals("Unknown root document element, beauti")) {
                infoLogger.severe("Error running file: " + fileName);
                infoLogger.severe(
                        "The file you just tried to run in BEAST is actually a BEAUti document.\n" +
                                "Although this uses XML, it is not a format that BEAST understands.\n" +
                                "These files are used by BEAUti to save and load your settings so that\n" +
                                "you can go back and alter them. To generate a BEAST file you must\n" +
                                "select the 'Generate BEAST File' option, either from the File menu or\n" +
                                "the button at the bottom right of the window.");

            } else {
                infoLogger.severe(
                        "\n******\nError in input BEAST XML file, " + fileName + ":\n" +
                        pxe.getMessage() +
                        "\n******\n"
                );
            }
            throw new RuntimeException("Terminate");
        } catch (RuntimeException rex) {

            if (rex.getMessage() != null && rex.getMessage().startsWith("The initial posterior is zero")) {
                infoLogger.severe("Error running file: " + fileName);
                infoLogger.severe(
                        "The initial model is invalid because state has a zero probability.\n\n" +
                                "If the log likelihood of the tree is -Inf, this may be because the\n" +
                                "initial, random tree is so large that it has an extremely bad\n" +
                                "likelihood which is being rounded to zero.\n\n" +
                                "Alternatively, it may be that the product of starting mutation rate\n" +
                                "and tree height is extremely small or extremely large. \n\n" +
                                "Finally, it may be that the initial state is incompatible with\n" +
                                "one or more 'hard' constraints (on monophyly or bounds on parameter\n" +
                                "values. This will result in Priors with zero probability.\n\n" +
                                "The individual components of the posterior are as follows:\n" +
                                rex.getMessage() + "\n" +
                                "For more information go to <http://beast.community>.");
            } else {
                // This call never returns as another RuntimeException exception is raised by
                // the error log handler???
                infoLogger.warning("Error running file: " + fileName);
                System.err.println("Fatal exception: " + rex.getMessage());
                rex.printStackTrace(System.err);
                System.err.flush();
            }
            throw new RuntimeException("Terminate");
        } catch (Exception ex) {
            infoLogger.warning("Error running file: " + fileName);
            infoLogger.severe("Fatal exception: " + ex.getMessage());
            System.err.println("Fatal exception: " + ex.getMessage());
            ex.printStackTrace(System.err);
            System.err.flush();
            throw new RuntimeException("Terminate");
        }

        try {
            dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate.releaseAllBeagleInstances();
            dr.evomodel.coalescent.basta.BeagleBastaLikelihoodDelegate.releaseAllBeagleBastaInstances();
        } catch (Throwable e) {
           throw new RuntimeException("Terminate");
        }
    }

    public static void centreLine(String line, int pageWidth) {
        int n = pageWidth - line.length();
        int n1 = n / 2;
        for (int i = 0; i < n1; i++) {
            System.out.print(" ");
        }
        System.out.println(line);
    }

    public static void printTitle() {
        System.out.println();
        centreLine("BEAST X " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("Bayesian Evolutionary Analysis Sampling Trees", 60);
        for (String creditLine : version.getCredits()) {
            centreLine(creditLine, 60);
        }
        System.out.println();

    }

    public static void printUsage(Arguments arguments) {

        arguments.printUsage("beast", "[<input-file-name>]");
        System.out.println();
        System.out.println("  Example: beast test.xml");
        System.out.println("  Example: beast -window test.xml");
        System.out.println("  Example: beast -help");
        System.out.println();
    }

    private static long updateSeedByRank(long seed, int rank) {
        return seed + 1000 * 1000 * rank;
    }

    //Main method
    public static void main(String[] args) throws java.io.IOException {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        Arguments arguments = new Arguments(
                new Arguments.Option[]{

                        new Arguments.Option("verbose", "Give verbose XML parsing messages"),
                        new Arguments.Option("warnings", "Show warning messages about BEAST XML file"),
                        new Arguments.Option("strict", "Fail on non-conforming BEAST XML file"),
                        new Arguments.Option("window", "Provide a console window"),
                        new Arguments.Option("options", "Display an options dialog"),
                        new Arguments.Option("working", "Change working directory to input file's directory"),
                        new Arguments.LongOption("seed", "Specify a random number generator seed"),
                        new Arguments.StringOption("prefix", "PREFIX", "Specify a prefix for all output log filenames"),
                        new Arguments.Option("overwrite", "Allow overwriting of log files"),
                        new Arguments.IntegerOption("errors", "Specify maximum number of numerical errors before stopping"),
                        new Arguments.IntegerOption("threads", "The maximum number of computational threads to use (default auto)"),
                        new Arguments.Option("fail_threads", "Exit with error on uncaught exception in thread"),
                        new Arguments.Option("ignore_versions", "Ignore mismatches between XML and BEAST versions"),
                        new Arguments.Option("java", "Use Java only, no native implementations"),
                        new Arguments.LongOption("tests", "The number of full evaluation tests to perform (default 1000)"),
                        new Arguments.RealOption("threshold", 0.0, Double.MAX_VALUE, "Full evaluation test threshold (default 0.1)"),
                        new Arguments.Option(OperatorSchedule.SHOW_OPERATORS, "Print transition kernel performance to file"),

                        new Arguments.Option("adaptation_off", "Don't adapt operator sizes"),
                        new Arguments.RealOption("adaptation_target", 0.0, 1.0, "Target acceptance rate for adaptive operators (default 0.234)"),

                        new Arguments.StringOption("pattern_compression", new String[]{"off", "unique", "ambiguous_constant", "ambiguous_all"},
                                false, "Site pattern compression mode - unique | ambiguous_constant | ambiguous_all (default unique)"),
                        new Arguments.RealOption("ambiguous_threshold", 0.0, 1.0, "Maximum proportion of ambiguous characters to allow compression (default 0.25)"),

                        new Arguments.Option("beagle", "Use BEAGLE library if available (default on)"),
                        new Arguments.Option("beagle_info", "BEAGLE: show information on available resources"),
                        new Arguments.Option("beagle_auto", "BEAGLE: automatically select fastest resource for analysis"),
                        new Arguments.StringOption("beagle_order", "order", "BEAGLE: set order of resource use"),
                        new Arguments.IntegerOption("beagle_instances", "BEAGLE: divide site patterns amongst instances"),
                        new Arguments.StringOption("beagle_multipartition", new String[]{"auto", "on", "off"},
                                false, "BEAGLE: use multipartition extensions if available (default auto)"),
                        new Arguments.Option("beagle_CPU", "BEAGLE: use CPU instance"),
                        new Arguments.Option("beagle_GPU", "BEAGLE: use GPU instance if available"),
                        new Arguments.Option("beagle_SSE", "BEAGLE: use SSE extensions if available"),
                        new Arguments.Option("beagle_SSE_off", "BEAGLE: turn off use of SSE extensions"),
                        new Arguments.Option("beagle_threading_off", "BEAGLE: turn off multi-threading for a CPU instance"),
                        new Arguments.IntegerOption("beagle_threads", 0, Integer.MAX_VALUE, "BEAGLE: manually set number of threads per CPU instance (default auto)"),
                        new Arguments.Option("beagle_cuda", "BEAGLE: use CUDA parallization if available"),
                        new Arguments.Option("beagle_opencl", "BEAGLE: use OpenCL parallization if available"),
                        new Arguments.Option("beagle_single", "BEAGLE: use single precision if available"),
                        new Arguments.Option("beagle_double", "BEAGLE: use double precision if available"),
                        new Arguments.Option("beagle_async", "BEAGLE: use asynchronous kernels if available"),
                        new Arguments.Option("beagle_low_memory", "BEAGLE: use lower memory pre-order traversal kernels"),
                        new Arguments.StringOption("beagle_extra_buffer_count", "buffer_count", "BEAGLE: reserve extra transition matrix buffers for convolutions"),
                        new Arguments.StringOption("beagle_scaling", new String[]{"default", "dynamic", "delayed", "always", "none"},
                                false, "BEAGLE: specify scaling scheme to use"),
                        new Arguments.Option("beagle_delay_scaling_off", "BEAGLE: don't wait until underflow for scaling option"),
                        new Arguments.LongOption("beagle_rescale", "BEAGLE: frequency of rescaling (dynamic scaling only)"),
                        new Arguments.Option("mpi", "Use MPI rank to label output"),

                        new Arguments.StringOption("particles", "FOLDER", "Specify a folder of particle start states"),

                        new Arguments.IntegerOption("mc3_chains", 1, Integer.MAX_VALUE, "number of chains"),
                        new Arguments.RealOption("mc3_delta", 0.0, Double.MAX_VALUE, "temperature increment parameter"),
                        new Arguments.RealArrayOption("mc3_temperatures", -1, "a comma-separated list of the hot chain temperatures"),
                        new Arguments.IntegerOption("mc3_swap", 1, Integer.MAX_VALUE, "frequency at which chains temperatures will be swapped"),
                        new Arguments.StringOption("mc3_scheme", "NAME", "Specify parallel tempering swap scheme"),

                        new Arguments.StringOption("load_state", "FILENAME", "Specify a filename to load a saved state from"),
                        new Arguments.StringOption("save_stem", "FILENAME", "Specify a stem for the filenames to save states to"),
                        new Arguments.LongOption("save_at", "Specify a state at which to save a state file"),
                        new Arguments.StringOption("save_time", "HH:mm:ss", "Specify a length of time after which to save a state file"),
                        new Arguments.LongOption("save_every", "Specify a frequency to save the state file"),
                        new Arguments.StringOption("save_state", "FILENAME", "Specify a filename to save state to"),
                        new Arguments.Option("full_checkpoint_precision", "Use hex-encoded doubles in checkpoint files"),
                        new Arguments.Option("force_resume", "Force resuming from a saved state"),

                        new Arguments.StringOption("citations_file", "FILENAME", "Specify a filename to write a citation list to"),
                        new Arguments.Option("citations_off", "Turn off writing citations to file"),
                        new Arguments.StringOption("plugins_dir", "FILENAME", "Specify a directory to load plugins from, multiple can be separated with ':' "),

                        new Arguments.Option("version", "Print the version and credits and stop"),
                        new Arguments.Option("help", "Print this information and stop"),
                });

        int argumentCount = 0;

        StringBuilder commandLine = new StringBuilder();
        for (int i = 0; i < args.length; i++) {
            if (i > 0) {
                commandLine.append(" ");
            }
            commandLine.append(args[i]);
        }
        System.setProperty("command_line", commandLine.toString());

        try {
            argumentCount = arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            System.out.println();
            System.out.println(ae.getMessage());
            System.out.println();
            printUsage(arguments);
            System.exit(1);
        }

        if (arguments.hasOption("version")) {
            printTitle();
        }

        if (arguments.hasOption("help")) {
            printUsage(arguments);
        }

        if (arguments.hasOption("version") || arguments.hasOption("help")) {
            System.exit(0);
        }

        List<String> additionalParsers = new ArrayList<String>();

        final boolean verbose = arguments.hasOption("verbose");
        if (verbose) {
            System.setProperty("verbose_output", Boolean.toString(true));
        }

        final boolean warnings = arguments.hasOption("warnings"); // if dev, then auto turn on, otherwise default to turn off
        if (warnings) {
            System.setProperty("show_warnings", Boolean.toString(true));
        }

        if (arguments.hasOption("ignore_versions")) {
            System.setProperty("ignore.versions", Boolean.toString(true));
        }

        final boolean strictXML = arguments.hasOption("strict");
        final boolean working = arguments.hasOption("working");
        String fileNamePrefix = null;
        boolean allowOverwrite = arguments.hasOption("overwrite");
        boolean useMPI = arguments.hasOption("mpi");

        long seed = MathUtils.getSeed();
        boolean useJava = false;

        int threadCount = -1;

        boolean usingSMC = false;

        boolean usingMC3 = false;
        double[] chainTemperatures = null;
        int swapChainsEvery = DEFAULT_SWAP_CHAIN_EVERY;
        MCMCMCOptions.SwapScheme swapScheme = MCMCMCOptions.SwapScheme.ORIGINAL_FLAVOR;

        if (arguments.hasOption("particles")) {
            System.setProperty("smc.particle_folder", arguments.getStringOption("particles"));
            usingSMC = true;

            System.setProperty("mcmc.evaluation.count", Long.toString(0));
        }

        // smc option is always run without GUI.
        final boolean window = !usingSMC && arguments.hasOption("window");
        final boolean options = !usingSMC && (arguments.hasOption("options") || (argumentCount == 0));

        if (!usingSMC) {
            if (arguments.hasOption("tests")) {
                long fullEvaluationCount = arguments.getLongOption("tests");
                System.setProperty("mcmc.evaluation.count", Long.toString(fullEvaluationCount));
            }

            if (arguments.hasOption("threshold")) {
                double fullEvaluationThreshold = arguments.getRealOption("threshold");
                System.setProperty("mcmc.evaluation.threshold", Double.toString(fullEvaluationThreshold));
            }

            if (arguments.hasOption("java")) {
                useJava = true;
            }

            if (arguments.hasOption("prefix")) {
                fileNamePrefix = arguments.getStringOption("prefix");
            }

            // ============= Evaluation approximations =============

            if (arguments.hasOption("pattern_compression")) {
                System.setProperty("patterns.compression", arguments.getStringOption("pattern_compression").toLowerCase());
            }
            if (arguments.hasOption("ambiguous_threshold")) {
                System.setProperty("patterns.threshold", String.valueOf(arguments.getRealOption("ambiguous_threshold")));
            }

            // ============= MC^3 settings =============

            int chainCount = 1;
            if (arguments.hasOption("mc3_chains")) {
                chainCount = arguments.getIntegerOption("mc3_chains");
            } else if (arguments.hasOption("mc3_temperatures")) {
                chainCount = 1 + arguments.getRealArrayOption("mc3_temperatures").length;
            }

            double delta = DEFAULT_DELTA;
            if (arguments.hasOption("mc3_delta")) {
                if (arguments.hasOption("mc3_temperatures")) {
                    System.err.println("Either the -mc3_delta or the -mc3_temperatures option should be used, not both");
                    System.err.println();
                    printUsage(arguments);
                    System.exit(1);
                }
                delta = arguments.getRealOption("mc3_delta");
            }

            chainTemperatures = new double[chainCount];
            chainTemperatures[0] = 1.0;
            if (arguments.hasOption("mc3_temperatures")) {
                double[] hotChainTemperatures = arguments.getRealArrayOption("mc3_temperatures");
                assert hotChainTemperatures.length == chainCount - 1;

                System.arraycopy(hotChainTemperatures, 0, chainTemperatures, 1, chainCount - 1);
            } else {
                for (int i = 1; i < chainCount; i++) {
                    chainTemperatures[i] = 1.0 / (1.0 + (delta * i));
                }
            }

            if (arguments.hasOption("mc3_swap")) {
                swapChainsEvery = arguments.getIntegerOption("mc3_swap");
            }

            if (arguments.hasOption("mc3_scheme")) {
                swapScheme = MCMCMCOptions.SwapScheme.parse(arguments.getStringOption("mc3_scheme"));
            }

            usingMC3 = chainCount > 1;
        }

        // ============= BEAGLE settings =============
        long beagleFlags = 0;

        boolean beagleShowInfo = arguments.hasOption("beagle_info");

        // if any beagle flag is specified then use beagle...
        //final boolean useBeagle = true;

        if (arguments.hasOption("beagle_auto")) {
            System.setProperty("beagle.resource.auto", Boolean.TRUE.toString());
        }
        if (arguments.hasOption("beagle_CPU")) {
            beagleFlags |= BeagleFlag.PROCESSOR_CPU.getMask();
        }
        if (arguments.hasOption("beagle_GPU")) {
            beagleFlags |= BeagleFlag.PROCESSOR_GPU.getMask();
        }
        if (arguments.hasOption("beagle_cuda")) {
            beagleFlags |= BeagleFlag.FRAMEWORK_CUDA.getMask();
        }
        if (arguments.hasOption("beagle_opencl")) {
            beagleFlags |= BeagleFlag.FRAMEWORK_OPENCL.getMask();
        }
        if (!arguments.hasOption("beagle_SSE_off")) {
            beagleFlags |= BeagleFlag.VECTOR_SSE.getMask();
        } else {
            beagleFlags |= BeagleFlag.VECTOR_NONE.getMask();
        }
        if (arguments.hasOption("beagle_thread_count")) {
            System.setProperty("beagle.thread.count", Integer.toString(arguments.getIntegerOption("beagle_thread_count")));
        }
        if (arguments.hasOption("beagle_threads")) {
            System.setProperty("beagle.thread.count", Integer.toString(arguments.getIntegerOption("beagle_threads")));
        }
        if (arguments.hasOption("beagle_threading_off")) {
            System.setProperty("beagle.thread.count", Integer.toString(1));
        }
//        if (arguments.hasOption("beagle_double")) {
//            beagleFlags |= BeagleFlag.PRECISION_DOUBLE.getMask();
//        }
        if (arguments.hasOption("beagle_single")) {
            beagleFlags |= BeagleFlag.PRECISION_SINGLE.getMask();
        } else {
            beagleFlags |= BeagleFlag.PRECISION_DOUBLE.getMask();
        }
        if (arguments.hasOption("beagle_async")) {
            beagleFlags |= BeagleFlag.COMPUTATION_ASYNCH.getMask();
        }

        if (arguments.hasOption("beagle_low_memory")) {
            beagleFlags |= BeagleFlag.PREORDER_TRANSPOSE_LOW_MEMORY.getMask();
        }

        if (arguments.hasOption("beagle_order")) {
            System.setProperty("beagle.resource.order", arguments.getStringOption("beagle_order"));
        }

        if (arguments.hasOption("beagle_instances")) {
            System.setProperty("beagle.instance.count", Integer.toString(arguments.getIntegerOption("beagle_instances")));
        }

        if (arguments.hasOption("beagle_scaling")) {
            System.setProperty("beagle.scaling", arguments.getStringOption("beagle_scaling"));
        }

        if (arguments.hasOption("beagle_delay_scaling_off")) {
            System.setProperty("beagle.delay.scaling", Boolean.FALSE.toString());
        }

        if (arguments.hasOption("beagle_rescale")) {
            System.setProperty("beagle.rescale", Long.toString(arguments.getLongOption("beagle_rescale")));
        }

        if (arguments.hasOption("beagle_multipartition")) {
            // attempt to force the use/non-use of the multipartition extensions, otherwise auto
            if (arguments.getStringOption("beagle_multipartition").toLowerCase().equals("on")) {
                System.setProperty("beagle.multipartition.extensions", Boolean.TRUE.toString());
            } else if (arguments.getStringOption("beagle_multipartition").toLowerCase().equals("off")) {
                System.setProperty("beagle.multipartition.extensions", Boolean.FALSE.toString());
            }
        }

        if (arguments.hasOption("beagle_extra_buffer_count")) {
            System.setProperty("beagle.extra.buffer.count", arguments.getStringOption("beagle_extra_buffer_count"));
        }

        // ============= Other settings =============
        if (arguments.hasOption("threads")) {
            // threadCount defaults to -1 unless the user specifies an option
            threadCount = arguments.getIntegerOption("threads");
            if (threadCount < 0) {
                printTitle();
                System.err.println("The the number of threads should be >= 0");
                System.exit(1);
            }
        }

        if (arguments.hasOption("fail_threads")) {

            Thread.UncaughtExceptionHandler handler = new Thread.UncaughtExceptionHandler() {
                @Override
                public void uncaughtException(Thread t, Throwable e) {
                    System.err.println("Error in thread " + t.getName() + ": ");
                    e.printStackTrace();
                    System.exit(-1);
                }
            };

            Thread.setDefaultUncaughtExceptionHandler(handler);
        }

        if (arguments.hasOption("seed")) {
            seed = arguments.getLongOption("seed");
            System.setProperty(BeastCheckpointer.CHECKPOINT_SEED, Long.toString(seed));
            if (seed <= 0) {
                printTitle();
                System.err.println("The random number seed should be > 0");
                System.exit(1);
            }
        }

        if (arguments.hasOption("adaptation_off")) {
            System.setProperty("mcmc.use_adaptation", Boolean.FALSE.toString());
        }
        if (arguments.hasOption("adaptation_target")) {
            System.setProperty("adaptation_target",
                    Double.toString(arguments.getRealOption("mcmc.adaptation_target")));
        }
        if (arguments.hasOption("plugins_dir")) {
            System.setProperty("beast.plugins.dir",
                     arguments.getStringOption("plugins_dir")+":"+System.getProperty("beast.plugins.dir"));
        }

        if (arguments.hasOption("citations_off")) {
            System.setProperty("output_citations", Boolean.FALSE.toString());
        } else {
            System.setProperty("output_citations", Boolean.TRUE.toString());
        }

        if (arguments.hasOption(OperatorSchedule.SHOW_OPERATORS)) {
            System.setProperty(OperatorSchedule.SHOW_OPERATORS, Boolean.TRUE.toString());
        }

        if (!usingSMC) {
            // ignore these other options

            if (arguments.hasOption("load_state")) {
                String stateFile = arguments.getStringOption("load_state");
                System.setProperty(BeastCheckpointer.LOAD_STATE_FILE, stateFile);
            }

            if (arguments.hasOption("save_at")) {
                long saveAt = arguments.getLongOption("save_at");
                System.setProperty(BeastCheckpointer.SAVE_STATE_AT, Long.toString(saveAt));
            }

            if (arguments.hasOption("save_time")) {
                String saveTime = arguments.getStringOption("save_time");
                System.setProperty(BeastCheckpointer.SAVE_STATE_TIME, saveTime);
            }

            if (arguments.hasOption("save_every")) {
                long saveEvery = arguments.getLongOption("save_every");
                System.setProperty(BeastCheckpointer.SAVE_STATE_EVERY, Long.toString(saveEvery));
            }

            if (arguments.hasOption("save_state")) {
                String stateFile = arguments.getStringOption("save_state");
                System.setProperty(BeastCheckpointer.SAVE_STATE_FILE, stateFile);
            }

            if (arguments.hasOption("save_stem")) {
                String stemName = arguments.getStringOption("save_stem");
                System.setProperty(BeastCheckpointer.SAVE_STEM, stemName);
            }

            if (arguments.hasOption("full_checkpoint_precision")) {
                System.setProperty(BeastCheckpointer.FULL_CHECKPOINT_PRECISION, "true");
            }

            if (arguments.hasOption("force_resume")) {
                System.setProperty("force.resume", Boolean.TRUE.toString());
            }

            if (arguments.hasOption("citations_file")) {
                String citationsFile = arguments.getStringOption("citations_file");
                System.setProperty("citations.filename", citationsFile);
            }

            if (useMPI) {
                String[] nullArgs = new String[0];
                try {
                    BeastMPI.Init(nullArgs);
                } catch (Exception e) {
                    throw new RuntimeException("Unable to access MPI.");
                }
                int rank = BeastMPI.COMM_WORLD.Rank();
                System.setProperty("mpi.rank.postfix", String.valueOf(rank));

            }

            String rankProp = System.getProperty("mpi.rank.postfix");
            if (rankProp != null) {
                int rank = Integer.valueOf(rankProp);
                seed = updateSeedByRank(seed, rank);
            }
        }

        int maxErrorCount = 0;
        if (arguments.hasOption("errors")) {
            maxErrorCount = arguments.getIntegerOption("errors");
            if (maxErrorCount < 0) {
                maxErrorCount = 0;
            }
        }

        BeastConsoleApp consoleApp = null;

        // don't include the pre-release commit in the window title bar
        String nameString = "BEAST X v" + version.getVersion();

        if (window) {
            System.setProperty("com.apple.macos.useScreenMenuBar", "true");
            System.setProperty("apple.laf.useScreenMenuBar", "true");
            System.setProperty("apple.awt.showGrowBox", "true");

            javax.swing.Icon icon = IconUtils.resize(IconUtils.getIcon(BeastMain.class, "images/beast.png"), 128, 128);

            String titleString = "<html>" +
                    "<div style=\"font: HelveticaNeue, Helvetica, Arial, sans-serif\">" +
                    "<p style=\"font-weight: 100; font-size: 42px\">BEAST X</p>" +
                    "<p style=\"font-weight: 200; font-size: 12px\">Bayesian Evolutionary Analysis Sampling Trees</p>" +
                    "<p style=\"font-weight: 300; font-size: 11px\">Version " + version.getVersionString() + ", " + version.getDateString() + "</p>" +
                    "</div></html>";

            String aboutString = "<html>" +
                    "<div style=\"font-family:HelveticaNeue-Light, 'Helvetica Neue Light', Helvetica, Arial, 'Lucida Grande',sans-serif; " +
                    "font-weight: 100;" +
                    "background: #35484F; color: #CBB944\">" +
                    "<center>" +
                    version.getHTMLCredits() +
                    "</div></center></div></html>";

            consoleApp = new BeastConsoleApp(nameString, titleString, aboutString,  icon);
            consoleApp.initialize();

        }

        printTitle();

        File inputFile = null;

        if (options && !beagleShowInfo) {

            String titleString = "<html>" +
                    "<div style=\"font-family: HelveticaNeue-Light, Helvetica, Arial, sans-serif\">" +
                    "<div style=\"font-family: HelveticaNeue-Thin; font-weight: 80; font-size: 42px\">BEAST X</div>" +
                    "<div style=\"font-weight: 200; font-size: 11px\">Bayesian Evolutionary Analysis Sampling Trees</div>" +
                    "<div style=\"font-weight: 300; font-size: 10px\">Version " + version.getVersionString() + ", " + version.getDateString() + "</div>" +
                    "<div style=\"font-weight: 300; font-size: 10px\"><a href=\"" + version.getBuildString() + "\">" +
                    version.getBuildString() + "</a></div>" +
                    "</div></html>";
            javax.swing.Icon icon = IconUtils.resize(IconUtils.getIcon(BeastMain.class, "images/beast.png"), 128, 128);

            BeastDialog dialog = new BeastDialog(new JFrame(), titleString, icon);

            dialog.setAllowOverwrite(allowOverwrite);
            dialog.setSeed(seed);

            //dialog.setUseBeagle(useBeagle);

            if (BeagleFlag.PROCESSOR_GPU.isSet(beagleFlags)) {
                dialog.setPreferBeagleGPU();
            }

            dialog.setPreferBeagleSSE(BeagleFlag.VECTOR_SSE.isSet(beagleFlags));

            if (BeagleFlag.PRECISION_SINGLE.isSet(beagleFlags)) {
                dialog.setPreferBeagleSingle();
            }

            if (!dialog.showDialog(nameString)) {
                return;
            }

            if (dialog.allowOverwrite()) {
                allowOverwrite = true;
            }

            seed = dialog.getSeed();
            threadCount = dialog.getThreadPoolSize();

            beagleShowInfo = dialog.showBeagleInfo();
            if (dialog.preferBeagleCPU()) {
                beagleFlags |= BeagleFlag.PROCESSOR_CPU.getMask();
            }
            if (dialog.preferBeagleSSE()) {
                beagleFlags |= BeagleFlag.VECTOR_SSE.getMask();
            } else {
                beagleFlags &= ~BeagleFlag.VECTOR_SSE.getMask();
            }
            if (dialog.preferBeagleGPU()) {
                beagleFlags |= BeagleFlag.PROCESSOR_GPU.getMask();
            }
            if (dialog.preferBeagleDouble()) {
                beagleFlags |= BeagleFlag.PRECISION_DOUBLE.getMask();
            }
            if (dialog.preferBeagleSingle()) {
                beagleFlags |= BeagleFlag.PRECISION_SINGLE.getMask();
            }
            System.setProperty("beagle.scaling", dialog.scalingScheme());

            inputFile = dialog.getInputFile();
            if (!beagleShowInfo && inputFile == null) {
                System.err.println("No input file specified");
                return;
            }

        }

        BeagleInfo.printVersionInformation();

        int[] versionNumbers = BeagleInfo.getVersionNumbers();
        if (versionNumbers.length != 0 && versionNumbers[0] < 4) {
            System.err.println("BEAST v" + BeastVersion.INSTANCE.getVersion() + " requires BEAGLE v4.0 or later.\n" +
                    "Please install or upgrade to the latest BEAGLE from https://beagle-dev.github.io/");
            throw new RuntimeException("Terminate");
        }

        if (beagleShowInfo) {
            BeagleInfo.printResourceList();
            return;
        }

        if (inputFile == null) {

            String[] args2 = arguments.getLeftoverArguments();

            if (args2.length > 1) {
                System.err.println("Unknown option: " + args2[1]);
                System.err.println();
                printUsage(arguments);
                return;
            }

            String inputFileName = null;

            if (args2.length > 0) {
                inputFileName = args2[0];
                inputFile = new File(inputFileName);
            }

            if (inputFileName == null) {
                // No input file name was given so throw up a dialog box...
                inputFile = Utils.getLoadFile("BEAST X " + version.getVersionString() + " - Select XML input file");
            }
        }

        if (inputFile != null && inputFile.getParent() != null && working) {
            System.setProperty("user.dir", inputFile.getParent());
        }

        if (window) {
            if (inputFile == null) {
                consoleApp.setTitle("null");
            } else {
                consoleApp.setTitle(inputFile.getName());
            }
        }

        if (useJava) {
            System.setProperty("java.only", "true");
        }

        if (fileNamePrefix != null && fileNamePrefix.trim().length() > 0) {
            System.setProperty("file.name.prefix", fileNamePrefix.trim());
        }

        if (allowOverwrite) {
            System.setProperty("log.allow.overwrite", "true");
        }

        additionalParsers.add("beagle");

        if (beagleFlags != 0) {
            System.setProperty("beagle.preferred.flags", Long.toString(beagleFlags));

        }

        if (threadCount >= 0) {
            System.setProperty("thread.count", String.valueOf(threadCount));
        }

        MathUtils.setSeed(seed);

        System.out.println("Random number seed: " + seed);

        try {
            new BeastMain(inputFile, consoleApp, maxErrorCount, verbose, warnings, strictXML, additionalParsers,
                    usingMC3 ? new MCMCMCOptions(chainTemperatures, swapChainsEvery, swapScheme) : null);
//                    usingMC3, chainTemperatures, swapChainsEvery);
        } catch (RuntimeException rte) {
            // The stack trace here is not useful
//            rte.printStackTrace(System.err);
            if (window) {
                System.out.println();
                System.out.println("BEAST has terminated with an error. Please select QUIT from the menu.");
                // logger.severe will throw a RTE but we want to keep the console visible
            } else {
                System.out.flush();
                System.err.flush();
                System.exit(1);
            }
        }

        if (useMPI) {
            BeastMPI.Finalize();
        }

        if (!usingMC3 && !window) { // DM: with MC3 the main thread gets here and terminates all threads prematurely 
            System.exit(0);
        }
    }
}


