package keser_master;

import Objects.Constants;
import org.biojava.nbio.core.sequence.DNASequence;

import java.text.DecimalFormat;
import java.util.List;

public class ToolMethods {
    public static DecimalFormat df_short = new DecimalFormat("0.00");
    public static DecimalFormat df = new DecimalFormat("0.0000");
    public static DecimalFormat df_long = new DecimalFormat("0.0000000");

    //calculates a value which determines how often a muatation results in a stop codon for single nucleotide mutations
    public static double[] getWeightedCountOfStopCodons_SNP() {
        double[] sum = new double[3];
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                for (int c = 0; c < 4; c++) {
                    if (Constants.isStopCodon(a, b, c)) continue;
                    for (int x = 0; x < 4; x++) {
                        for (int pos = 0; pos < 3; pos++) {
                            String NewCodon = "";
                            int basePrev;
                            if (pos == 0) {
                                NewCodon = Constants.Bases[x] + Constants.Bases[b] + Constants.Bases[c];
                                basePrev = a;
                            } else if (pos == 1) {
                                NewCodon = Constants.Bases[a] + Constants.Bases[x] + Constants.Bases[c];
                                basePrev = b;
                            } else {
                                NewCodon = Constants.Bases[a] + Constants.Bases[b] + Constants.Bases[x];
                                basePrev = c;
                            }
                            if (!Constants.isStopCodon(NewCodon)) continue;
                            double count = 1;
                            if (MainClass.nucleotideApriori != null) {
                                count = count * MainClass.nucleotideApriori.getValue(basePrev);
                            }
                            if (MainClass.tripletApriori != null) {
                                count = count * MainClass.tripletApriori.getValue(a, b, c);
                            }
                            sum[pos] = sum[pos] + count;
                        }
                    }
                }
            }
        }
        return sum;
    }

    //calculates a value which determines how often a muatation results in a stop codon for frameshift muatations
    public static double[] getWeightedCountOfStopCodons_Shift() {
        double[] sum = new double[2];
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                for (int c = 0; c < 4; c++) {
                    if (Constants.isStopCodon(a, b, c)) continue;
                    for (int x = 0; x < 4; x++) {
                        for (int dir = 0; dir < 2; dir++) {
                            String NewCodon;
                            if (dir == 0) {
                                NewCodon = Constants.Bases[x] + Constants.Bases[a] + Constants.Bases[b];
                            } else {
                                NewCodon = Constants.Bases[b] + Constants.Bases[c] + Constants.Bases[x];
                            }
                            if (!Constants.isStopCodon(NewCodon)) continue;
                            double count = 1;
                            if (MainClass.nucleotideTransition != null) {
                                if (dir == 0) {
                                    count = count * MainClass.nucleotideTransition.getValue(x, a);
                                } else {
                                    count = count * MainClass.nucleotideTransition.getValue(c, x);
                                }
                            }
                            if (MainClass.tripletTransition != null) {
                                if (dir == 0) {
                                    count = count * MainClass.tripletTransition.getValueLeft(a, b, c, x);
                                } else {
                                    count = count * MainClass.tripletTransition.getValueRight(a, b, c, x);
                                }

                            }
                            sum[dir] = sum[dir] + count;
                        }
                    }
                }
            }
        }
        return sum;
    }

    public static double getWeightedStopCodonFrequency_Overall() {
        double[] SNP = getWeightedCountOfStopCodons_SNP();
        double[] Shift = getWeightedCountOfStopCodons_Shift();
        return SNP[0] + SNP[1] + SNP[2] + Shift[0] + Shift[1];
    }

    //calculates the ratio of coding tripletts to stop-tripletts
    public static double getStopCodonFrequency(List<DNASequence> sequences, boolean inReadingFrame, int offset) {
        int triplettcount = 0;
        int stopcodoncount = 0;
        for (DNASequence seq : sequences) {
            String sequence = seq.getSequenceAsString();
            if (offset != 0) {
                sequence = sequence.substring(offset);
            }
            int stepsize = 1;
            if (inReadingFrame) {
                stepsize = 3;
            }
            //count stopcodons
            for (int i = 0; i < sequence.length() - 2; i = i + stepsize) {
                String codon = sequence.substring(i, i + 3);
                triplettcount++;
                if (Constants.isStopCodon(codon)) {
                    stopcodoncount++;
                }
            }

        }
        return (double) stopcodoncount / (double) triplettcount;
    }

    //Elementwise Difference between 2 4x4 matrices
    public static double[][] MatrixDiff(double[][] M1, double[][] M2) {
        double[][] Erg = new double[4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                Erg[i][j] = M1[i][j] - M2[i][j];
            }
        }
        double sum = 0;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                sum += Math.abs(Erg[i][j]);
            }
        }
        System.out.println("Gesamtdiferenz: " + df.format(sum));
        PrintMatrix(Erg);
        return Erg;
    }

    //Prints the 4x4 transition matrix in the console
    public static void PrintMatrix(double[][] M) {
        System.out.println("Vertikal: s(n) horizontal: s(n+1)");
        System.out.println("--- T ------- C ------- A ------- G");
        System.out.println("T " + df.format(M[0][0]) + " -- " + df.format(M[1][0]) + " -- " + df.format(M[2][0]) + " -- " + df.format(M[3][0]));
        System.out.println("C " + df.format(M[0][1]) + " -- " + df.format(M[1][1]) + " -- " + df.format(M[2][1]) + " -- " + df.format(M[3][1]));
        System.out.println("A " + df.format(M[0][2]) + " -- " + df.format(M[1][2]) + " -- " + df.format(M[2][2]) + " -- " + df.format(M[3][2]));
        System.out.println("G " + df.format(M[0][3]) + " -- " + df.format(M[1][3]) + " -- " + df.format(M[2][3]) + " -- " + df.format(M[3][3]));
    }

    //Prints a 4x4x4 matrix in the console
    public static void PrintMatrix(double[][][] M) {
        String[] N = Constants.Bases;
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                System.out.println(N[a] + N[b] + N[0] + ": " + df.format(M[a][b][0]) + " " + N[a] + N[b] + N[1] + ": " + df.format(M[a][b][1]) + " " + N[a] + N[b] + N[2] + ": " + df.format(M[a][b][2]) + " " + N[a] + N[b] + N[3] + ": " + df.format(M[a][b][3]) + " ");
            }
        }
    }

    //Prints NA Weightings
    public static void PrintMatrix(double[] M) {
        System.out.println("T: " + df.format(M[0])+ "A-priori: "+ df.format(M[0]*0.25));
        System.out.println("C: " + df.format(M[1])+ "A-priori: "+ df.format(M[1]*0.25));
        System.out.println("A: " + df.format(M[2])+ "A-priori: "+ df.format(M[2]*0.25));
        System.out.println("G: " + df.format(M[3])+ "A-priori: "+ df.format(M[3]*0.25));
    }

    //Prints a 4x4x4 matrix in the console
    public static void PrintTripletTable(double[][][] M, boolean shortDecimal) {
        DecimalFormat d = shortDecimal ? df_short : df;
        String[] N = Constants.Bases;
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                System.out.println(N[a] + N[b] + N[0] + " & " + d.format(M[a][b][0]) + " & " + N[a] + N[b] + N[1] + " & " + d.format(M[a][b][1]) + " & " + N[a] + N[b] + N[2] + " & " + d.format(M[a][b][2]) + " & " + N[a] + N[b] + N[3] + " & " + d.format(M[a][b][3]));
            }
        }
    }

    public static double[][][] calculateZScores(double[][][] matrix) {
        double[][][] zScores = new double[4][4][4];
        //mean value
        double mean = 0;
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                for (int c = 0; c < 4; c++) {
                    mean += matrix[a][b][c];
                }
            }
        }
        mean = mean/64;
        //standard deviation
        double sigma = 0;
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                for (int c = 0; c < 4; c++) {
                    sigma = sigma + Math.pow(matrix[a][b][c]-mean, 2);
                }
            }
        }
        sigma = Math.sqrt(sigma / 63); //divide my N-1

        //z-scores
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                for (int c = 0; c < 4; c++) {
                    zScores[a][b][c] = (matrix[a][b][c] - mean) / sigma;
                }
            }
        }
        System.out.println("Mean: "+df.format(mean)+" Sigma: "+df.format(sigma));
        return zScores;
    }

    public static double[][][][][][] calculateZScores(double[][][][][][] matrix) {
        double[][][][][][] zScores = new double[4][4][4][4][4][4];
        //mean value
        double mean = 0;
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                for (int c = 0; c < 4; c++) {
                    for (int x = 0; x < 4; x++) {
                        for (int y = 0; y < 4; y++) {
                            for (int z = 0; z < 4; z++) {
                                mean += matrix[a][b][c][x][y][z];
                            }
                        }
                    }
                }
            }
        }
        mean = mean/(64*64);
        //standard deviation
        double sigma = 0;
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                for (int c = 0; c < 4; c++) {
                    for (int x = 0; x < 4; x++) {
                        for (int y = 0; y < 4; y++) {
                            for (int z = 0; z < 4; z++) {
                                sigma = sigma + Math.pow(matrix[a][b][c][x][y][z]-mean, 2);
                            }
                        }
                    }
                }
            }
        }
        sigma = Math.sqrt(sigma / (64*64 - 1)); //divide my N-1

        //z-scores
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                for (int c = 0; c < 4; c++) {
                    for (int x = 0; x < 4; x++) {
                        for (int y = 0; y < 4; y++) {
                            for (int z = 0; z < 4; z++) {
                                zScores[a][b][c][x][y][z] = (matrix[a][b][c][x][y][z] - mean) / sigma;
                            }
                        }
                    }
                }
            }
        }
        System.out.println("Mean: "+df.format(mean)+" Sigma: "+df.format(sigma));
        return zScores;
    }

    //Prints a 4x4x4 matrix in the console
    public static void printTT2Table(double[][][][][][] M) {
        DecimalFormat d =  df_short;
        String[] N = Constants.Bases;
        String firstCol = "";
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                for (int c = 0; c < 4; c++) {
                    firstCol = "\\hline\\multirow{4}{*}{" + N[a] + N[b] + N[c] +"}";
                    for (int x = 0; x < 4; x++) {
                       System.out.println(firstCol +"&" + N[x]+"&"+
                                       d.format(M[a][b][c][x][0][0])  + "&" +
                                       d.format(M[a][b][c][x][0][1])  + "&" +
                                       d.format(M[a][b][c][x][0][2])  + "&" +
                                       d.format(M[a][b][c][x][0][3])  + "&" +

                                       d.format(M[a][b][c][x][1][0])  + "&" +
                                       d.format(M[a][b][c][x][1][1])  + "&" +
                                       d.format(M[a][b][c][x][1][2])  + "&" +
                                       d.format(M[a][b][c][x][1][3])  + "&" +

                                       d.format(M[a][b][c][x][2][0])  + "&" +
                                       d.format(M[a][b][c][x][2][1])  + "&" +
                                       d.format(M[a][b][c][x][2][2])  + "&" +
                                       d.format(M[a][b][c][x][2][3])  + "&" +

                                       d.format(M[a][b][c][x][3][0])  + "&" +
                                       d.format(M[a][b][c][x][3][1])  + "&" +
                                       d.format(M[a][b][c][x][3][2])  + "&" +
                                       d.format(M[a][b][c][x][3][3])  + " \\\\");
                       firstCol = "";
                    }
                }
            }
        }
    }

}
