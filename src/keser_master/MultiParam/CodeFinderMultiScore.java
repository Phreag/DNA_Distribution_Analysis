package keser_master.MultiParam;

import keser_master.MultiParam.code.IGreedyCode;
import keser_master.MultiParam.code.MultiScoreCode;
import keser_master.Objects.GeneCode;
import keser_master.StabilityCalculator;
import keser_master.ToolMethods;

import java.io.FileWriter;
import java.io.IOException;

public class CodeFinderMultiScore extends CodeFinderAbstract{

    void calculateNaturalCodeValues(){
        MultiScoreCode naturalCode = new MultiScoreCode(GeneCode.naturalCode);
        GeneCode natural = new GeneCode();
        StabilityCalculator stCalc=new StabilityCalculator(natural);

        naturalCode.setMS1(stCalc.get_BaseDeviation(1));
        naturalCode.setMS2(stCalc.get_BaseDeviation(2));
        naturalCode.setMS3(stCalc.get_BaseDeviation(3));

        naturalCode.setrMS(stCalc.get_ShiftDeviation(1));
        naturalCode.setlMS(stCalc.get_ShiftDeviation(2));

        this.naturalCode = naturalCode;
    }

    @Override
    Thread getThreadedCalculator() {
        return new Thread(() -> {
            GeneCode g=new GeneCode();
            StabilityCalculator S=new StabilityCalculator(g);
            while (true){
                int currentCode=getNextValue();
                if(currentCode>=codesWithoutValues.size())break;
                String[] rCode=codesWithoutValues.get(currentCode);
                MultiScoreCode codeData = new MultiScoreCode(rCode);
                g.changeCode(rCode);

                codeData.setMS1(S.get_BaseDeviation(1));
                codeData.setMS2(S.get_BaseDeviation(2));
                codeData.setMS3(S.get_BaseDeviation(3));

                codeData.setrMS(S.get_ShiftDeviation(1));
                codeData.setlMS(S.get_ShiftDeviation(2));

                codeData.calculateGreedyScore((MultiScoreCode)naturalCode);

                ValueBuffer[currentCode] = codeData;
            }
        });
    }

    void writeValueLine(FileWriter writer, int Number, IGreedyCode c){
        //Writes Value Line To File
        MultiScoreCode code = (MultiScoreCode) c;
        try {
            writer.write(Number+" "+String.join("~", c.getCode())+ " Score: "+c.getGreedy_score() +"["
                    + ToolMethods.df_short.format(code.getMS1()) +"|"
                    + ToolMethods.df_short.format(code.getMS2()) +"|"
                    + ToolMethods.df_short.format(code.getMS3()) +"|"
                    + ToolMethods.df_short.format(code.getrMS()) +"|"
                    + ToolMethods.df_short.format(code.getlMS()) +"]\n");
        } catch (IOException e) {
            System.out.println("FileWriter Error");
        }
    }
}
