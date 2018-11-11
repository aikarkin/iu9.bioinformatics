class AlignmentConfiguration {
    private char compound;
    private String firstSeqFile, secondSeqFile;
    private int open = -10;
    private int extend = -1;
    private String alignmentFile;

    public void setAlignmentFile(String alignmentFile) {
        this.alignmentFile = alignmentFile;
    }


    public AlignmentConfiguration(char compound, String firstSeqFile, String secondSeqFile, int open, int extend) {
        this.compound = compound;
        this.firstSeqFile = firstSeqFile;
        this.secondSeqFile = secondSeqFile;
//        this.alignmentFile = alignmentFile;
        this.open = open;
        this.extend = extend;
    }

    public AlignmentConfiguration(char compound, String firstSeqFile, String secondSeqFile, String alignmentFile) {
        this.compound = compound;
        this.firstSeqFile = firstSeqFile;
        this.secondSeqFile = secondSeqFile;
        this.alignmentFile = alignmentFile;
    }


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

    public int getOpen() {
        return open;
    }


}