import org.apache.commons.cli.*;

import javax.naming.ConfigurationException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Optional;
import java.util.function.BiFunction;

import static java.lang.Math.min;

public class RunSequenceAlignment {
    private static final int MAX_CHARS_PER_LINE = 80;
    private static final int OPEN = -10, EXTEND = -1;

    @FunctionalInterface
    interface ScoreFunction {
        int apply(int [][]arr, int i, int j);
    }

    public static void main(String[] args) {
        try {
            String firstSeq, secondSeq;

            // Обрабатываем аргумменты командной строки ( помощью библиотеки Apache Commons CLI)
            CommandLine cmd = initCommandLine(args);

            if(cmd != null) {
                // По распознаным опциям формируем конфигурацию будущего выравнивания
                AlignmentConfiguration conf = alignConfigFromCmd(cmd);

                // Инициализируем поток вывода
                PrintWriter outWriter = conf.getAlignmentFile() == null
                        ? new PrintWriter(System.out, true)
                        : new PrintWriter(conf.getAlignmentFile());

                firstSeq = readInputSequence(conf.getFirstSeqFile());
                secondSeq = readInputSequence(conf.getSecondSeqFile());

                BiFunction<Character, Character, Integer> scoreFunction;

                // если нужно выравнить последовательность аминокислот, используем скоринг функцию blosum62, если нуклеотидов - dnaFul
                // иначе - выдаем ошибку.
                if (conf.getCompound() == 'a' && NWUtils.isAminoAcidsSequence(firstSeq) && NWUtils.isAminoAcidsSequence(secondSeq)) {
                    System.out.println("Type: amino acid");
                    scoreFunction = NWUtils::blosum62;
                } else if (conf.getCompound() == 'n' && NWUtils.isNucleotideSequence(firstSeq) && NWUtils.isNucleotideSequence(secondSeq)) {
                    System.out.println("Type: nucleotide");
                    scoreFunction = NWUtils::dnaFull;
                } else {
                    System.err.println("[error] Invalid input sequence");
                    return;
                }

                // запускам выравнивание с заданными параметрами
                alignSequences(firstSeq, secondSeq, OPEN, EXTEND, outWriter, scoreFunction);
            }

        } catch (ConfigurationException e) {
            String msg = e.getMessage();
            System.err.printf("[error] Invalid alignment configuration %s.\n", (msg == null || msg.length() == 0) ? "" : ": " + msg);
        } catch (IOException e) {
            String msg = e.getMessage();
            e.printStackTrace();
            System.err.printf("[error] Unable to read file%s.\n", (msg == null || msg.length() == 0) ? "" : ": " + msg);
        }
    }

    private static void initMatrices(int[][] matrix_m, int[][] matrix_i, int[][] matrix_d, int n, int m, int open, int extend) {
        int inf = 2 * open + (n + m) * extend + 1;

        matrix_m[0][0] = 0;
        matrix_d[0][0] = matrix_i[0][0] = inf;

        for (int i = 1; i < n; i++) {
            matrix_m[i][0] = matrix_d[i][0] = inf;
            matrix_i[i][0] = open + (i - 1) * extend;
        }

        for (int j = 1; j < m; j++) {
            matrix_m[0][j] = matrix_d[0][j] = inf;
            matrix_i[0][j] = open + (j - 1) * extend;
        }
    }

    private static int max(int a, int b, int c) {
        return Math.max(a, Math.max(b, c));
    }

    public static void alignSequences(String firstSeq, String secondSeq, int open, int extend, PrintWriter out, BiFunction<Character, Character, Integer> scoreFunction) {
        int n = firstSeq.length() + 1,
            m = secondSeq.length() + 1;

        int[][] matrix_m = new int[n][m],
                matrix_i = new int[n][m],
                matrix_d = new int[n][m];

        initMatrices(matrix_m, matrix_i, matrix_d, n, m, open, extend);

        for (int i = 1; i < n; i++) {
            for (int j = 1; j < m; j++) {
                int score_m = scoreFunction.apply(firstSeq.charAt(i - 1), secondSeq.charAt(j - 1));
                matrix_m[i][j] = max(matrix_m[i - 1][j - 1] + score_m, matrix_i[i - 1][j - 1] + score_m, matrix_d[i - 1][j - 1] + score_m);
                matrix_i[i][j] = max(matrix_i[i][j - 1] + extend, matrix_m[i][j - 1] + open, matrix_d[i][j - 1] + open);
                matrix_d[i][j] = max(matrix_i[i - 1][j] + extend, matrix_m[i - 1][j] + open, matrix_d[i - 1][j] + open);
            }
        }


        StringBuilder firstBuilder = new StringBuilder();
        StringBuilder secondBuilder = new StringBuilder();
        int i = n - 1, j = m - 1;

        while(i > 0 || j > 0) {
            int optimalScore = max(matrix_i[i][j], matrix_d[i][j], matrix_m[i][j]);

            // match/mismatch
            if (i > 0 && j > 0 && optimalScore == matrix_m[i][j]) {
                System.out.print("match");
                firstBuilder.insert(0, firstSeq.charAt(i - 1));
                secondBuilder.insert(0, secondSeq.charAt(j - 1));
                i--; j--;
            // delete
            } else if(i > 0 && optimalScore == matrix_d[i][j]) {
                System.out.print("delete");
                firstBuilder.insert(0, firstSeq.charAt(i - 1));
                secondBuilder.insert(0, '_');
                i--;
            // insert
            } else if(j > 0 && optimalScore == matrix_i[i][j]) {
                System.out.print("insert");
                firstBuilder.insert(0, '_');
                secondBuilder.insert(0, secondSeq.charAt(j - 1));
                j--;
            }
            System.out.println(" " + optimalScore);
        }

        int score = max(matrix_i[n - 1][m - 1], matrix_d[n - 1][m - 1], matrix_m[n - 1][m - 1]);

        printAlignmentAndScore(score, firstBuilder.toString(), secondBuilder.toString(), out);
    }

    private static void printTable(PrintWriter out, int[][] score) {
        int m = score[0].length;
        int n = score.length;

        out.println("Score matrix: \n");

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                out.printf(" %s ", score[i][j]);
            }
            out.println();
        }
    }


    private static void printAlignmentAndScore(int score, String firstSeq, String secondSeq, PrintWriter out) {
        int n = firstSeq.length(),
                k = (n % MAX_CHARS_PER_LINE == 0) ? n / MAX_CHARS_PER_LINE : n / MAX_CHARS_PER_LINE + 1;

        String[] sequences = {firstSeq, secondSeq};

        out.println("Score: " + score);
        out.println();

        for (int i = 0; i < k; i++) {
            for(int si = 0; si < 2; si++) {
                out.printf("seq%d: ", si + 1);
                for (int j = MAX_CHARS_PER_LINE * i; j < min(n, MAX_CHARS_PER_LINE * (i + 1)); j++) {
                    out.print(sequences[si].charAt(j));
                }
                out.println();
            }
            out.println();
        }

    }

    private static String readInputSequence(String filepath) throws IOException {
        return new String(Files.readAllBytes(Paths.get(filepath)));
    }

    private static CommandLine initCommandLine(String[] args) {
        Options cmdOptions = new Options();

        cmdOptions.addOption(
            Option.builder("h")
                .longOpt("help")
                .desc("Print help message")
                .hasArg(false)
                .build()
        );
        cmdOptions.addOption(
                Option.builder("i")
                        .longOpt("input")
                        .desc("Two input files with sequences, that should be aligned")
                        .hasArgs()
                        .numberOfArgs(2)
                        .required()
                        .build()
        );
        cmdOptions.addOption(
                Option.builder("c")
                        .longOpt("compound")
                        .desc("Type of organic compound. Available values: 'a' - amino acids, 'n' - nucleotide  Default value: 'a'.")
                        .type(String.class)
                        .hasArg()
                        .build()
        );
        cmdOptions.addOption(
                Option.builder("o")
                        .longOpt("open")
                        .desc("Fine of open gap.")
                        .hasArg()
                        .numberOfArgs(1)
                        .type(Integer.class)
                        .required()
                        .build()
        );
        cmdOptions.addOption(
                Option.builder("e")
                        .longOpt("extend")
                        .desc("Fine of gap extend.")
                        .hasArg()
                        .numberOfArgs(1)
                        .type(Integer.class)
                        .required()
                        .build()
        );
        cmdOptions.addOption(
                Option.builder()
                        .argName("a")
                        .desc("If provided match output alignment file path with result alignment and score.")
                        .longOpt("output")
                        .hasArg()
                        .type(String.class)
                        .build()
        );

        CommandLine cmd = null;
        try {
             cmd = new DefaultParser().parse(cmdOptions, args);

        } catch (ParseException e) {
            String msg = e.getMessage();
            System.err.printf("[error] Invalid command line arguments%s.\n", (msg == null || msg.length() == 0) ? "" : ": " + msg);
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp( "align", cmdOptions, true);
        }
        return cmd;
    }

    private static AlignmentConfiguration alignConfigFromCmd(CommandLine cmd) throws ConfigurationException {
        String compound, seq1File, seq2File;
        int open, extend;

        String[] seqFiles = cmd.getOptionValues('i');

        if(seqFiles.length != 2) {
            throw new ConfigurationException("Invalid number of input files.");
        }

        seq1File = seqFiles[0];
        seq2File = seqFiles[1];

        compound = cmd.hasOption('c') ? cmd.getOptionValue('c') : "a";

        if(!compound.equals("a") && !compound.equals("n")) {
            throw new ConfigurationException("Invalid compound '" + compound + "'. Compound may has values 'a' or 'n'." );
        }

        open = parseInt(cmd.getOptionValue('o')).orElse(-10);
        extend = parseInt(cmd.getOptionValue("g")).orElse(-1);

        if(open > 0) {
            throw new ConfigurationException("Invalid open value. It should be negative integer");
        }

        if(extend > 0) {
            throw new ConfigurationException("Invalid extend value. It should be negative integer");
        }

        AlignmentConfiguration conf = new AlignmentConfiguration(compound.charAt(0), seq1File, seq2File, open, extend);

        if(cmd.hasOption('a')) {
            conf.setAlignmentFile(cmd.getOptionValue('a'));
        }

        return conf;
    }

    private static Optional<Integer> parseInt(String toParse) {
        try {
            return Optional.of(Integer.parseInt(toParse));
        } catch (NumberFormatException e) {
            return Optional.empty();
        }
    }
}
