import org.apache.commons.cli.*;

import javax.naming.ConfigurationException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Optional;
import java.util.function.BiFunction;

import static java.lang.Math.max;
import static java.lang.Math.min;

public class RunSequenceAlignment {
    private static final int MAX_CHARS_PER_LINE = 80;

    public static void main(String[] args) {
        try {
            String firstSeq, secondSeq;

            CommandLine cmd = initCommandLine(args);
            AlignmentConfiguration conf = alignConfigFromCmd(cmd);

            PrintWriter outWriter = conf.getAlignmentFile() == null
                    ? new PrintWriter(System.out, true)
                    : new PrintWriter(conf.getAlignmentFile());

            firstSeq = readInputSequence(conf.getFirstSeqFile());
            secondSeq = readInputSequence(conf.getSecondSeqFile());

            BiFunction<Character, Character, Integer> scoreFunction;

            if(conf.getCompound() == 'a' && NWUtils.isAminoAcidsSequence(firstSeq) && NWUtils.isAminoAcidsSequence(secondSeq)) {
                System.out.println("Type: amino acid");
                scoreFunction = NWUtils::blosum62;
            } else if(conf.getCompound() == 'n' && NWUtils.isNucleotideSequence(firstSeq) && NWUtils.isNucleotideSequence(secondSeq)) {
                System.out.println("Type: nucleotide");
                scoreFunction = NWUtils::dnaFull;
            } else {
                System.err.println("[error] Invalid input sequence");
                return;
            }

            alignSequences(firstSeq, secondSeq, conf.getGap(), outWriter, scoreFunction);

        } catch (ParseException | ConfigurationException e) {
            String msg = e.getMessage();
            System.err.printf("[error] Invalid command line arguments%s.\n", (msg == null || msg.length() == 0) ? "" : ": " + msg);
        } catch (IOException e) {
            String msg = e.getMessage();
            e.printStackTrace();
            System.err.printf("[error] Unable to read file%s.\n", (msg == null || msg.length() == 0) ? "" : ": " + msg);
        }
    }

    private static void alignSequences(String firstSeq, String secondSeq, int gap, PrintWriter out, BiFunction<Character, Character, Integer> scoreFunction) {
        int match, delete, insert, n = firstSeq.length(), m = secondSeq.length();
        int[][] score = new int[n][m];

        // Этап 1. Формирование матрицы скора для каждой пары префиксов A, B
        score[0][0] = 0;

        for (int i = 0; i < n; i++) {
            score[i][0] = gap * i;
        }

        for (int j = 0; j < m; j++) {
            score[0][j] = gap * j;
        }

        for (int i = 1; i < n; i++) {
            for (int j = 1; j < m; j++) {
                // скор совпадения - диагональный элемент + скор функция для текущей клетки
                match = score[i - 1][j - 1] + scoreFunction.apply(firstSeq.charAt(i), secondSeq.charAt(j));
                // скор удаления - левый элемент + штраф за гэп
                delete = score[i - 1][j] + gap;
                // скор вставки - верхний элемент + штраф за гэп
                insert = score[i][j - 1] + gap;
                score[i][j] = max(max(match, delete), insert);
            }
        }

        // Этап 2. Восстанавливаем по матрице выравнивание
//        StringBuilder firstBuilder = new StringBuilder(), secondBuilder = new StringBuilder();
        StringBuilder firstBuilder = new StringBuilder();
        StringBuilder secondBuilder = new StringBuilder();
        int i = n - 1, j = m - 1;

        while(i > 0 || j > 0) {
            if(i > 0 && j > 0 && score[i][j] == score[i - 1][j - 1] + scoreFunction.apply(firstSeq.charAt(i), secondSeq.charAt(j))) {
                firstBuilder.insert(0, firstSeq.charAt(i));
                secondBuilder.insert(0, secondSeq.charAt(j));
                i--; j--;
            } else if(i > 0 && score[i][j] == score[i - 1][j] + gap) {
                firstBuilder.insert(0, firstSeq.charAt(i));
                secondBuilder.insert(0, '_');
                i--;
            } else if(j > 0 && score[i][j] == score[i][j - 1] + gap) {
                firstBuilder.insert(0, '_');
                secondBuilder.insert(0, secondSeq.charAt(j));
                j--;
            }
        }

        System.out.println("Length of aligned sequences: " + firstBuilder.toString().length());

        printAlignmentAndScore(score[n - 1][m - 1], firstBuilder.toString(), secondBuilder.toString(), out);
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

    private static CommandLine initCommandLine(String[] args) throws ParseException {
        Options cmdOptions = new Options();

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
                Option.builder("g")
                        .longOpt("gap")
                        .desc("Fine of gap in miss-matches.")
                        .hasArg()
                        .numberOfArgs(1)
                        .type(Integer.class)
                        .required()
                        .build()
        );
        cmdOptions.addOption(
                Option.builder()
                        .argName("o")
                        .longOpt("output")
                        .hasArg()
                        .type(String.class)
                        .build()
        );

        return new DefaultParser().parse(cmdOptions, args);
    }

    private static AlignmentConfiguration alignConfigFromCmd(CommandLine cmd) throws ConfigurationException {
        String compound, seq1File, seq2File;
        int gap;

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

        gap = parseInt(cmd.getOptionValue('g')).orElse(1);

        if(gap > 0) {
            throw new ConfigurationException("Invalid gap value. It should be negative integer");
        }

        AlignmentConfiguration conf = new AlignmentConfiguration(compound.charAt(0), seq1File, seq2File, gap);

        if(cmd.hasOption('o')) {
            conf.setAlignmentFile(cmd.getOptionValue('o'));
        }

        return conf;
    }

    private static class AlignmentConfiguration {
        private char compound;
        private String firstSeqFile, secondSeqFile;

        public void setAlignmentFile(String alignmentFile) {
            this.alignmentFile = alignmentFile;
        }

        private String alignmentFile;

        public AlignmentConfiguration(char compound, String firstSeqFile, String secondSeqFile, int gap) {
            this.compound = compound;
            this.firstSeqFile = firstSeqFile;
            this.secondSeqFile = secondSeqFile;
            this.alignmentFile = alignmentFile;
            this.gap = gap;
        }

        public AlignmentConfiguration(char compound, String firstSeqFile, String secondSeqFile, String alignmentFile) {
            this.compound = compound;
            this.firstSeqFile = firstSeqFile;
            this.secondSeqFile = secondSeqFile;
            this.alignmentFile = alignmentFile;
        }

        private int gap = -10;

        public char getCompound() {
            return compound;
        }

        public String getFirstSeqFile() {
            return firstSeqFile;
        }

        public String getSecondSeqFile() {
            return secondSeqFile;
        }

        public String getAlignmentFile() {
            return alignmentFile;
        }

        public int getGap() {
            return gap;
        }
    }

    private static Optional<Integer> parseInt(String toParse) {
        try {
            return Optional.of(Integer.parseInt(toParse));
        } catch (NumberFormatException e) {
            return Optional.empty();
        }
    }
}
