class AlignmentConfiguration {
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