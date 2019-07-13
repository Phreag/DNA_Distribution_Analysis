package keser_master.Objects;

import keser_master.MainClass;

public class WeightLoop {
    private NucleotideApriori na;
    private TripletApriori ta;
    private NucleotideTransition nt;
    private TripletTransition tt;
    private int currentpos = 0;
    private boolean excludeNT = false;

    public WeightLoop(NucleotideApriori na, TripletApriori ta, TripletTransition tt) {
        this(na,ta,null,tt);
        excludeNT=true;
    }


    public WeightLoop(NucleotideApriori na, TripletApriori ta, NucleotideTransition nt, TripletTransition tt) {
        this.na = na;
        this.ta = ta;
        this.nt = nt;
        this.tt = tt;
    }

    private void setWeightsInMainClass(int pos){
        boolean[] wts = {false, false, false, false};
        if(pos>=8){
            pos=pos-8;
            wts[0]=true;
        }
        if(pos>=4){
            pos=pos-4;
            wts[1]=true;
        }
        if(pos>=2){
            pos=pos-2;
            wts[2]=true;
        }
        if (pos>=1){
            wts[3]=true;
        }
        MainClass.setWeightings(wts[2] ? na : null, wts[1] ? ta : null, wts[0] ? nt : null, wts[3] ? tt : null);
    }

    //returns false if finished.
    public boolean moveNext(){
        if(excludeNT){
            if (currentpos > 7) return false;
        }else{
            if (currentpos > 15) return false;
        }
        setWeightsInMainClass(currentpos);
        currentpos++;
        return true;
    }
}
