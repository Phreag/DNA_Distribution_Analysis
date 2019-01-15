package Objects;


public class TripletTransition {
    private double [][][][][][] tripletTransition;
    private double [][][][][] tripletTransitionNucleotide;

    public TripletTransition(double[][][][][][] tripletTransition, double[][][][][] tripletTransitionNucleotide) {
        this.tripletTransition = tripletTransition;
        this.tripletTransitionNucleotide = tripletTransitionNucleotide;
    }

    public double getValue(int a, int b, int c, int x, int y, int z){
        return tripletTransition[a][b][c][x][y][z];
    }

    //old: Front
    public double getValueLeft(int a, int b, int c, int x){
        return tripletTransitionNucleotide[a][b][c][x][0];
    }

    //old: After
    public double getValueRight(int a, int b, int c, int x){
        return tripletTransitionNucleotide[a][b][c][x][1];
    }

    public double[][][][][] getDataNucleotide(){
        return tripletTransitionNucleotide;
    }

    public double[][][][][][] getDataTriplet(){
        return tripletTransition;
    }



}
