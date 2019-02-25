package Objects;


public class TripletTransition {
    private double [][][][][][] tripletTransition;
    private double [][][][][][] tripletTransitionClean;
    private double [][][][][] tripletTransitionNucleotide;


    public TripletTransition(double[][][][][][] tripletTransition, double[][][][][] tripletTransitionNucleotide) {
        this.tripletTransition = tripletTransition;
        this.tripletTransitionNucleotide = tripletTransitionNucleotide;
    }

    public double getValue(int a, int b, int c, int x, int y, int z){
        if(tripletTransitionClean == null){
            return tripletTransition[a][b][c][x][y][z];
        }
        return tripletTransitionClean[a][b][c][x][y][z];
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

    public void cleanTriplettApriori(TripletApriori TA){
        tripletTransitionClean = new double[4][4][4][4][4][4];
        //normalize for each of the first tripletts
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                for (int c = 0; c < 4; c++) {
                    double sum = 0;
                    //divide by TA
                    for (int x = 0; x < 4; x++) {
                        for (int y = 0; y < 4; y++) {
                            for (int z = 0; z < 4; z++) {
                                tripletTransitionClean[a][b][c][x][y][z] = tripletTransition[a][b][c][x][y][z]/TA.getValue(a,b,c);
                                sum += tripletTransitionClean[a][b][c][x][y][z];
                            }
                        }
                    }
                    sum = sum/64;
                    //normalize
                    for (int x = 0; x < 4; x++) {
                        for (int y = 0; y < 4; y++) {
                            for (int z = 0; z < 4; z++) {
                                tripletTransitionClean[a][b][c][x][y][z] = tripletTransitionClean[a][b][c][x][y][z]/sum;
                            }
                        }
                    }

                }
            }
        }

        //normalize for each of the second tripletts
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                for (int c = 0; c < 4; c++) {
                    double sum = 0;
                    //divide by TA
                    for (int x = 0; x < 4; x++) {
                        for (int y = 0; y < 4; y++) {
                            for (int z = 0; z < 4; z++) {
                                tripletTransitionClean[x][y][z][a][b][c] = tripletTransition[x][y][z][a][b][c]/TA.getValue(a,b,c);
                                sum += tripletTransitionClean[x][y][z][a][b][c];
                            }
                        }
                    }
                    sum = sum/64;
                    //normalize
                    for (int x = 0; x < 4; x++) {
                        for (int y = 0; y < 4; y++) {
                            for (int z = 0; z < 4; z++) {
                                tripletTransitionClean[x][y][z][a][b][c] = tripletTransitionClean[x][y][z][a][b][c]/sum;
                            }
                        }
                    }

                }
            }
        }

    }

}
