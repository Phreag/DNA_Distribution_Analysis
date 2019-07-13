package keser_master;

import keser_master.Objects.*;

import java.util.HashMap;
import java.util.Map;

public class MarkovChainForStopCodonsCalculator {
    private GeneCode code;
    private TripletApriori tripletApriori;
    private NucleotideTransition nucleotideTransition;
    private TripletTransition tripletTransition;
    private boolean useNucleotideTransition = true;
    private Map<String, Double> averageDistanceSums = new HashMap<>();
    private Map<String, Double> cachedValues = new HashMap<>();
    private int maxdepth = 2000;
    private boolean maxdepth_return_zero = false;
    private String searchCodon;

    private MarkovChainForStopCodonsCalculator(GeneCode code, TripletTransition tripletTransition, NucleotideTransition nucleotideTransition, TripletApriori tripletApriori,int maxdepth, boolean maxdepth_return_zero, String searchCodon){
        this.code = code;
        this.tripletTransition = tripletTransition;
        this.nucleotideTransition=nucleotideTransition;
        this.tripletApriori = tripletApriori;
        useNucleotideTransition= tripletTransition==null;
        this.maxdepth = maxdepth;
        this.maxdepth_return_zero = maxdepth_return_zero;
        this.searchCodon = searchCodon;
        calculateAverageDistSums();
    }

    public MarkovChainForStopCodonsCalculator(GeneCode code, TripletTransition tripletTransition, TripletApriori tripletApriori) {
        this(code,tripletTransition,null,tripletApriori,2000,false,null);
    }

    public MarkovChainForStopCodonsCalculator(GeneCode code, NucleotideTransition nucleotideTransition, TripletApriori tripletApriori) {
        this(code,null,nucleotideTransition,tripletApriori,2000,false,null);
    }
    public MarkovChainForStopCodonsCalculator(GeneCode code, TripletTransition tripletTransition, TripletApriori tripletApriori, int maxdepth, String searchCodon) {
        this(code,tripletTransition,null,tripletApriori,maxdepth,false,searchCodon);
    }

    public MarkovChainForStopCodonsCalculator(GeneCode code, NucleotideTransition nucleotideTransition, TripletApriori tripletApriori,int maxdepth, boolean maxdepth_return_zero) {
        this(code,null,nucleotideTransition,tripletApriori,maxdepth,maxdepth_return_zero,null);
    }

    public MarkovChainForStopCodonsCalculator(GeneCode code, TripletTransition tripletTransition, TripletApriori tripletApriori, int maxdepth, boolean maxdepth_return_zero){
        this(code,tripletTransition,null,tripletApriori,maxdepth,maxdepth_return_zero,null);
    }

    private boolean isSearchedCodon(String codon){
        if(searchCodon != null){
            return codon.equalsIgnoreCase(searchCodon);
        }else{
            return Constants.isStopCodon(codon);
        }
    }

    public double getChainLengthAfterMuatation() {
        double averageDistanceSum = 0;
        double count = 0;
        for (int i = 0; i < 4; i++) {
            String a = Constants.Bases[i];
            for (int j = 0; j < 4; j++) {
                String b = Constants.Bases[j];
                for (int k = 0; k < 4; k++) {
                    String c = Constants.Bases[k];
                    String codon1 = a + b + c;
                    String Amino = code.getAminoAcid(codon1);
                    if (Amino.length() != 3) continue; //Filters Stop Codons as Origin
                    for (int m = 0; m < 4; m++) {
                        String x = Constants.Bases[m];
                        String[] newCodons = new String[]{b + c + x, x + a + b};
                        for (String codon2 : newCodons) {
                            double weight = tripletApriori.getValue(i,j,k);
                            double averageDist = averageDistanceSums.get(codon2);
                            averageDistanceSum += (averageDist * weight);
                            count++;
                        }
                    }
                }
            }
        }
        double average = averageDistanceSum / count;
        //System.out.println("Average Distance to Stop Codon: " + ToolMethods.df.format(average));
        return average;
    }

    private void calculateAverageDistSums() {
        for (int i = 0; i < 4; i++) {
            String a = Constants.Bases[i];
            for (int j = 0; j < 4; j++) {
                String b = Constants.Bases[j];
                for (int k = 0; k < 4; k++) {
                    String c = Constants.Bases[k];
                    String codon = a + b + c;
                    averageDistanceSums.put(codon, getDistanceToStopCodonRec(0, i, j, k));
                    //System.out.println("Distance " + codon + ": " + ToolMethods.df.format(averageDistanceSums.get(codon)));
                }
            }
        }
    }

    private double getWeight(int a, int b, int c, int x, int y, int z){
        if(useNucleotideTransition){
            return nucleotideTransition.getValue(a,b)*
                    nucleotideTransition.getValue(b,c)*
                    nucleotideTransition.getValue(c,x)*
                    nucleotideTransition.getValue(x,y)*
                    nucleotideTransition.getValue(y,z);
        }else{
            return tripletTransition.getValue(a,b,c,x,y,z);
        }
    }


    private double getDistanceToStopCodonRec(int depth, int a, int b, int c) {
        if (depth >= maxdepth){
            if(maxdepth_return_zero){
                return 0;
            }
            return maxdepth;
        }
        String codon = Constants.Bases[a] + Constants.Bases[b] + Constants.Bases[c];
        if (cachedValues.containsKey(codon + "_" + depth)) {
            return cachedValues.get(codon + "_" + depth);
        }
        if (isSearchedCodon(codon)) return depth;
        double[][] distArray = new double[64][2];
        double weightsSum = 0;
        //0 = weight
        //1 = dist
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    int index = k + (4 * j) + (16 * i);
                    //double weight = 1;
                    double weight = getWeight(a,b,c,i,j,k);
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
}
