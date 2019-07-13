package keser_master.MultiParam.code;

public class MultiCharacteristicsCode implements IGreedyCode {
    //Holds a GeneCode and its calculated GMS Scores on different amino acid characteristics
    private final String[] code;
    private double GMS_Polar;
    private double GMS_Hydro;
    private double GMS_MolVol;
    private double GMS_Mr;
    private double GMS_pKa;
    private double GMS_pKb;
    private double GMS_PI;

    private double greedy_score;

    public MultiCharacteristicsCode(String[]  code){
        this.code = code;
    }

    public String[]  getCode() {
        return code;
    }

    public double getGMS_Polar() {
        return GMS_Polar;
    }

    public void setGMS_Polar(double GMS_Polar) {
        this.GMS_Polar = GMS_Polar;
    }

    public double getGMS_Hydro() {
        return GMS_Hydro;
    }

    public void setGMS_Hydro(double GMS_Hydro) {
        this.GMS_Hydro = GMS_Hydro;
    }

    public double getGMS_MolVol() {
        return GMS_MolVol;
    }

    public void setGMS_MolVol(double GMS_MolVol) {
        this.GMS_MolVol = GMS_MolVol;
    }

    public double getGMS_Mr() {
        return GMS_Mr;
    }

    public void setGMS_Mr(double GMS_Mr) {
        this.GMS_Mr = GMS_Mr;
    }

    public double getGMS_pKa() {
        return GMS_pKa;
    }

    public void setGMS_pKa(double GMS_pKa) {
        this.GMS_pKa = GMS_pKa;
    }

    public double getGMS_pKb() {
        return GMS_pKb;
    }

    public void setGMS_pKb(double GMS_pKb) {
        this.GMS_pKb = GMS_pKb;
    }

    public double getGMS_PI() {
        return GMS_PI;
    }

    public void setGMS_PI(double GMS_PI) {
        this.GMS_PI = GMS_PI;
    }

    public double getGreedy_score() {
        return greedy_score;
    }

    public void calculateGreedyScore(MultiCharacteristicsCode natCode) {

        if(natCode.getGMS_Polar() <= GMS_Polar ||
                natCode.getGMS_Hydro() <= GMS_Hydro ||
                natCode.getGMS_MolVol() <= GMS_MolVol ||
                natCode.getGMS_Mr() <= GMS_Mr ||
                natCode.getGMS_pKa() <= GMS_pKa ||
                natCode.getGMS_pKb() <= GMS_pKb ||
                natCode.getGMS_PI() <= GMS_PI) {
            //at least one value is greater than the value from the natural code
            greedy_score = Math.max(0, GMS_Polar - natCode.getGMS_Polar())
                    + Math.max(0, GMS_Hydro - natCode.getGMS_Hydro())
                    + Math.max(0, GMS_MolVol - natCode.getGMS_MolVol())
                    + Math.max(0, GMS_Mr - natCode.getGMS_Mr())
                    + Math.max(0, GMS_pKa - natCode.getGMS_pKa())
                    + Math.max(0, GMS_pKb - natCode.getGMS_pKb())
                    + Math.max(0, GMS_PI - natCode.getGMS_PI());
        }else{
           //all values are smaller than value from the natural code (negative score)
            greedy_score = GMS_Polar + GMS_Hydro + GMS_MolVol + GMS_Mr + GMS_pKa + GMS_pKb + GMS_PI
                    - natCode.getGMS_Polar()
                    - natCode.getGMS_Hydro()
                    - natCode.getGMS_MolVol()
                    - natCode.getGMS_Mr()
                    - natCode.getGMS_pKa()
                    - natCode.getGMS_pKb()
                    - natCode.getGMS_PI();
        }


    }
}
