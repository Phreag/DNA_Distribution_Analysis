package Objects;

public class SequenceStatsCalculator {
    private String Sequence;
    private final boolean readInReadingFrame;
    private int offset;
    //rawData contains the Counts
    private int[][][][][][] rawData = new int[4][4][4][4][4][4];

    /**
     * when readInReadingFrame = true the statistics will be calulated in steps by 3
     * the offset is only effective if 	readInReadingFrame = true. it is used to read an alternative reading frame
     */
    public SequenceStatsCalculator(boolean readInReadingFrame, int offset) {
        this.readInReadingFrame = readInReadingFrame;
        this.offset = offset;
    }

    public void processSequence(String Sequence) {
        this.Sequence = Sequence;
        UpdateRawData();
    }

    private void UpdateRawData() {
        int step = 1;
        if (readInReadingFrame) {
            step = 3;
        }
        if (offset != 0) {
            this.Sequence = this.Sequence.substring(offset);
        }

        for (int i = 0; i < Sequence.length() - 6; i = i + step) {
            int[] nucleotides = {-1, -1, -1, -1, -1, -1};
            int remaining = 6;
            for (int n = 0; n < 6; n++) {
                try {
                    nucleotides[n] = Constants.getNucleotideNumber(Sequence.charAt(i + n));
                    remaining--;
                } catch (IllegalArgumentException e) {
                    break;
                }
            }
            if (remaining == 0) {
                rawData[nucleotides[0]][nucleotides[1]][nucleotides[2]][nucleotides[3]][nucleotides[4]][nucleotides[5]]++;
            }
        }
    }

    public NucleotideApriori getNucleotide_aPriori() {
        int[] BaseCount = new int[4];
        int overallCount = 0;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    for (int l = 0; l < 4; l++) {
                        for (int m = 0; m < 4; m++) {
                            for (int n = 0; n < 4; n++) {
                                if (readInReadingFrame) {
                                    BaseCount[i] += rawData[i][j][k][l][m][n];
                                    BaseCount[j] += rawData[i][j][k][l][m][n];
                                    BaseCount[k] += rawData[i][j][k][l][m][n];
                                    overallCount += rawData[i][j][k][l][m][n] * 3;
                                } else {
                                    BaseCount[i] += rawData[i][j][k][l][m][n];
                                    overallCount += rawData[i][j][k][l][m][n];
                                }
                            }
                        }
                    }
                }
            }
        }

        double[] NA = new double[4];
        for (int i = 0; i < 4; i++) {
            NA[i] = ((double) BaseCount[i] / (double) overallCount) * 4;
        }
        System.out.println("Nucleotide_aPriori Average: " + (NA[0] + NA[1] + NA[2] + NA[3]) / 4);
        return new NucleotideApriori(NA);
    }


    public TripletApriori getTriplet_aPriori() {
        int[][][] tripletCount = new int[4][4][4];
        int overallCount = 0;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    for (int l = 0; l < 4; l++) {
                        for (int m = 0; m < 4; m++) {
                            for (int n = 0; n < 4; n++) {
                                //if(Constants.isStopCodon(i,j,k))continue;
                                tripletCount[i][j][k] += rawData[i][j][k][l][m][n];
                                overallCount += rawData[i][j][k][l][m][n];
                            }
                        }
                    }
                }
            }
        }
        double[][][] triplet_aPriori = new double[4][4][4];
        double sum = 0;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    triplet_aPriori[i][j][k] = ((double) tripletCount[i][j][k] * 64) / (double) overallCount;
                    sum += triplet_aPriori[i][j][k];
                }
            }
        }
        sum = sum / 64;
        System.out.println("Triplet_aPriori Average: " + sum);
        return new TripletApriori(triplet_aPriori);
    }

    public NucleotideTransition getNucleotideTransition() {
        int[][] nucleotideTransitionCount = new int[4][4];
        long overallCount = 0;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    for (int l = 0; l < 4; l++) {
                        for (int m = 0; m < 4; m++) {
                            for (int n = 0; n < 4; n++) {
                                if (readInReadingFrame) {
                                    //if(Constants.isStopCodon(i,j,k))continue;
                                    nucleotideTransitionCount[i][j] += rawData[i][j][k][l][m][n];
                                    nucleotideTransitionCount[j][k] += rawData[i][j][k][l][m][n];
                                    nucleotideTransitionCount[k][l] += rawData[i][j][k][l][m][n];
                                    overallCount += rawData[i][j][k][l][m][n]*3;
                                } else {
                                    nucleotideTransitionCount[i][j] += rawData[i][j][k][l][m][n];
                                    overallCount += rawData[i][j][k][l][m][n];
                                }

                            }
                        }
                    }
                }
            }
        }
        double[][] nucleotideTransition = new double[4][4];
        double sum = 0;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                nucleotideTransition[i][j] = ((double) nucleotideTransitionCount[i][j] / (double) overallCount) * 16;
                sum += nucleotideTransition[i][j];
            }
        }
        sum = sum / 16;
        System.out.println("BaseTransition Average: " + sum);
        return new NucleotideTransition(nucleotideTransition);
    }

    public TripletTransition gettripletTransition() {
        int overallCount = 0;
        int[][][][][] tripletTransitionNucleotideCount = new int[4][4][4][4][2];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    for (int l = 0; l < 4; l++) {
                        for (int m = 0; m < 4; m++) {
                            for (int n = 0; n < 4; n++) {
                                overallCount += rawData[i][j][k][l][m][n];
                                tripletTransitionNucleotideCount[i][j][k][l][1] += rawData[i][j][k][l][m][n];
                                tripletTransitionNucleotideCount[l][m][n][k][0] += rawData[i][j][k][l][m][n];
                            }
                        }
                    }
                }
            }
        }
        double[][][][][][] tripletTransition = new double[4][4][4][4][4][4];
        double[][][][][] tripletTransitionNucleotide = new double[4][4][4][4][2];
        double sumFront = 0;
        double sumAfter = 0;
        double sum = 0;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    for (int l = 0; l < 4; l++) {
                        for (int m = 0; m < 4; m++) {
                            for (int n = 0; n < 4; n++) {
                                tripletTransition[i][j][k][l][m][n] = (double) rawData[i][j][k][l][m][n] * 4096 / (double) overallCount;
                                sum += tripletTransition[i][j][k][l][m][n];
                                tripletTransitionNucleotide[i][j][k][l][1] += (double) rawData[i][j][k][l][m][n] * 256 / (double) overallCount;
                                tripletTransitionNucleotide[l][m][n][k][0] += (double) rawData[i][j][k][l][m][n] * 256 / (double) overallCount;
                                sumFront += tripletTransitionNucleotide[i][j][k][l][1];
                                sumAfter += tripletTransitionNucleotide[l][m][n][k][0];
                            }
                        }
                    }
                }
            }
        }
        sum = sum / 4096;
        sumFront=sumFront/256;
        sumAfter=sumAfter/256;
        System.out.println("Triplet Transition Average (Nucleotide): Front:"+sumFront+" After:"+sumAfter);
        System.out.println("Triplet Transition Average: " + sum);
        return new TripletTransition(tripletTransition, tripletTransitionNucleotide);
    }
}
