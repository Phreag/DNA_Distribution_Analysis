package Objects;

public class NucleotideApriori {
    private double[] nucleotide_aPriori;

    public NucleotideApriori(double[] nucleotide_aPriori) {
        this.nucleotide_aPriori = nucleotide_aPriori;
    }
    public double getValue(int nucleotide){
        return nucleotide_aPriori[nucleotide];
    }
}
