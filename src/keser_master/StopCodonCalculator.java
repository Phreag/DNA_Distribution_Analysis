package keser_master;

import Objects.Constants;
import Objects.GeneCode;

import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map;

public class StopCodonCalculator {
    private GeneCode code;
    private double[][][] tripletAprioriWeights;
    private double[][] baseTransitionWeights;
    static DecimalFormat df = new DecimalFormat("0.0000");
    private String[] Bases = {"T", "C", "A", "G"};
    private Map<String, Double> averageDistanceSums = new HashMap<>();
    private Map<String, Double> cachedValues = new HashMap<>();

    public StopCodonCalculator(GeneCode code, double[][] baseTransitionWeights, double[][][] tripletAprioriWeights) {
        this.code = code;
        this.baseTransitionWeights = baseTransitionWeights;
        this.tripletAprioriWeights = tripletAprioriWeights;
        calculateAverageDistSums();
    }


    public double getChainLengthAfterMuatation() {
        double averageDistanceSum = 0;
        double count = 0;

        for (int i = 0; i < 4; i++) {
            String a = Bases[i];
            for (int j = 0; j < 4; j++) {
                String b = Bases[j];
                for (int k = 0; k < 4; k++) {
                    String c = Bases[k];
                    String codon1 = a + b + c;
                    String Amino = code.getAminoAcid(codon1);
                    if (Amino.length() != 3) continue; //Filters Stop Codons as Origin
                    for (int m = 0; m < 4; m++) {
                        String x = Bases[m];
                        String[] newCodons = new String[]{b + c + x, x + a + b};
                        for (String codon2 : newCodons) {
                            double weight = tripletAprioriWeights[i][j][k];
                            double averageDist = averageDistanceSums.get(codon2);
                            averageDistanceSum += (averageDist * weight);
                            count++;
                        }
                    }
                }
            }
        }

        double average = averageDistanceSum / count;
        System.out.println("Average Distance to Stop Codon: " + df.format(average));
        return average;
    }

    private void calculateAverageDistSums() {
        for (int i = 0; i < 4; i++) {
            String a = Bases[i];
            for (int j = 0; j < 4; j++) {
                String b = Bases[j];
                for (int k = 0; k < 4; k++) {
                    String c = Bases[k];
                    String codon = a + b + c;
                    averageDistanceSums.put(codon, getDistanceToStopCodonRec(0, i, j, k));
                    System.out.println("Distance " + codon + ": " + df.format(averageDistanceSums.get(codon)));
                }
            }
        }
    }


    private double getDistanceToStopCodonRec(int depth, int a, int b, int c) {
        int maxdepth = 2000;
        if (depth >= maxdepth) return maxdepth;
        String codon = Bases[a] + Bases[b] + Bases[c];
        if (cachedValues.containsKey(codon + "_" + depth)) {
            return cachedValues.get(codon + "_" + depth);
        }

        String Amino = code.getAminoAcid(codon);
        if (Amino.length() != 3) return depth;
        double[][] distArray = new double[64][2];
        double weightsSum = 0;
        //0 = weight
        //1 = dist
        for (int i = 0; i < 4; i++) {
            String x = Bases[i];
            for (int j = 0; j < 4; j++) {
                String y = Bases[j];
                for (int k = 0; k < 4; k++) {
                    String z = Bases[k];
                    int index = k + (4 * j) + (16 * i);
                    //double weight = 1;
                    double weight = (baseTransitionWeights[c][i]) * (baseTransitionWeights[i][j]) * (baseTransitionWeights[j][k]);
                    distArray[index][0] = weight;
                    distArray[index][1] = getDistanceToStopCodonRec(depth + 1, i, j, k);
                    weightsSum += weight;
                }
            }
        }
        double averageWeight = (weightsSum / 64);
        double dist = 0;
        for (int i = 0; i < 64; i++) {
            dist = dist + (distArray[i][0] / averageWeight) * (distArray[i][1] / 64);
        }
        cachedValues.put(codon + "_" + depth, dist);
        return dist;
    }



    /*
    Für jedes Triplett
        Für jede Mutation dieses Tripletts
            Wie viele Basen kommen danach im Leseraster im Durchschitt bis ein Stoppcodon im Leseraster erscheint
     */
}
