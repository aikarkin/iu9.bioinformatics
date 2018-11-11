import org.junit.jupiter.api.Test;

import java.io.PrintWriter;
import java.util.function.BiFunction;

public class Tests {
    @Test
    public void test1() {
        int open = -10, extend = -1;
        BiFunction<Character, Character, Integer> scoreFunc = (a, b) -> (a == b) ? 5 : 4;
        String seq1 = "ACGT", seq2 = "ACGGCTT";
        RunSequenceAlignment.alignSequences(seq1, seq2, open, extend, new PrintWriter(System.out, true), scoreFunc);
    }
}
