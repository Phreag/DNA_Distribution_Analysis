package Objects;

import static Objects.Constants.*;

public class Constants_Object {
    private int mode = 0;
    /*
    0: polarReq
    1: hydropathy
    2: molVol
    3: molWeight
    4: pKa
    5: pKb
    6: pI
     */

    public void nextMode(){
        mode++;
        if(mode > 6){
            mode = 0;
        }
    }
    public void reset(){
        mode = 0;
    }


    public double getSqareDifference(String oldAminoAcid, String newAminoAcid) {
        double diff;
        switch (mode){
            case 0:
                double polar1 = getPolarReq(oldAminoAcid);
                double polar2 = getPolarReq(newAminoAcid);
                diff = polar1 - polar2;
                return Math.pow(diff, 2);
            case 1:
                double hydro1 = getHydropathy(oldAminoAcid);
                double hydro2 = getHydropathy(newAminoAcid);
                diff = hydro1 - hydro2;
                return Math.pow(diff, 2);
            case 2:
                double value1 = getMolVol(oldAminoAcid);
                double value2 = getMolVol(newAminoAcid);
                diff = value1 - value2;
                return Math.pow(diff, 2);
            case 3:
                value1 = getMolWeight(oldAminoAcid);
                value2 = getMolWeight(newAminoAcid);
                diff = value1 - value2;
                return Math.pow(diff, 2);
            case 4:
                value1 = getpKa(oldAminoAcid);
                value2 = getpKa(newAminoAcid);
                diff = value1 - value2;
                return Math.pow(diff, 2);
            case 5:
                value1 = getpKb(oldAminoAcid);
                value2 = getpKb(newAminoAcid);
                diff = value1 - value2;
                return Math.pow(diff, 2);
            case 6:
                value1 = getpI(oldAminoAcid);
                value2 = getpI(newAminoAcid);
                diff = value1 - value2;
                return Math.pow(diff, 2);
        }
        throw new IllegalStateException();
    }

}
