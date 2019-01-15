package Objects;

public class TripletApriori {
    private double[][][] triplet_aPriori;
    public TripletApriori(double[][][] triplet_aPriori){
        this.triplet_aPriori=triplet_aPriori;
    }

    //return the apriori weighting for the triplett abc
    public double getValue(int a, int b, int c){
        return triplet_aPriori[a][b][c];
    }

    public double[][][] getData(){
        return triplet_aPriori;
    }
}
