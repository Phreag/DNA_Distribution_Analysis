package keser_master.MultiParam.code;

public class MultiScoreCode implements IGreedyCode{
    private final String[] code;
    private double greedy_score;

    private double MS1;
    private double MS2;
    private double MS3;
    private double rMS;

    public double getMS1() {
        return MS1;
    }

    public void setMS1(double MS1) {
        this.MS1 = MS1;
    }

    public double getMS2() {
        return MS2;
    }

    public void setMS2(double MS2) {
        this.MS2 = MS2;
    }

    public double getMS3() {
        return MS3;
    }

    public void setMS3(double MS3) {
        this.MS3 = MS3;
    }

    public double getrMS() {
        return rMS;
    }

    public void setrMS(double rMS) {
        this.rMS = rMS;
    }

    public double getlMS() {
        return lMS;
    }

    public void setlMS(double lMS) {
        this.lMS = lMS;
    }

    private double lMS;

    public MultiScoreCode(String[] code) {
        this.code = code;
    }

    @Override
    public double getGreedy_score() {
        return greedy_score;
    }

    @Override
    public String[] getCode() {
        return code;
    }

    public void calculateGreedyScore(MultiScoreCode natCode) {
        if(natCode.getMS1() <= MS1 ||
                natCode.getMS2() <= MS2 ||
                natCode.getMS3() <= MS3 ||
                natCode.getrMS() <= rMS ||
                natCode.getlMS() <= lMS ) {
            //at least one value is greater than the value from the natural code
            greedy_score = Math.max(0, MS1 - natCode.getMS1())
                    + Math.max(0, MS2 - natCode.getMS2())
                    + Math.max(0, MS3 - natCode.getMS3())
                    + Math.max(0, rMS - natCode.getrMS())
                    + Math.max(0, lMS - natCode.getlMS());
        }else{
            //all values are smaller than value from the natural code (negative score)
            greedy_score = MS1 + MS2 + MS3 + rMS + lMS
                    - natCode.getMS1()
                    - natCode.getMS2()
                    - natCode.getMS3()
                    - natCode.getrMS()
                    - natCode.getlMS();
        }
    }
}
