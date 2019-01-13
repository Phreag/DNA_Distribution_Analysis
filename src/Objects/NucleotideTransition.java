package Objects;

public class NucleotideTransition {
    private double[][] nucleotideTransition;

    public NucleotideTransition(double[][] nucleotideTransition) {
        this.nucleotideTransition = nucleotideTransition;
    }

    public double getValue(int a, int b){
        return nucleotideTransition[a][b];
    }
}
