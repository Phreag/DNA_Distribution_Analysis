package keser_master.MultiParam;

import keser_master.MultiParam.code.IGreedyCode;
import keser_master.MultiParam.code.MultiCharacteristicsCode;
import keser_master.Objects.Constants_Object;
import keser_master.Objects.GeneCode;
import keser_master.StabilityCalculator;
import keser_master.ToolMethods;

import java.io.FileWriter;
import java.io.IOException;

public class CodeFinderMultiCharacteristics extends CodeFinderAbstract{

    void calculateNaturalCodeValues(){
        MultiCharacteristicsCode naturalCode = new MultiCharacteristicsCode(GeneCode.naturalCode);
        GeneCode natural = new GeneCode();
        Constants_Object constants = new Constants_Object();
        StabilityCalculator stCalc=new StabilityCalculator(natural, constants);
        constants.reset();
        //0: polarReq
        naturalCode.setGMS_Polar(stCalc.getGMS());
        //1: hydropathy
        constants.nextMode();
        naturalCode.setGMS_Hydro(stCalc.getGMS());
        //2: molVol
        constants.nextMode();
        naturalCode.setGMS_MolVol(stCalc.getGMS());
        //3: molWeight
        constants.nextMode();
        naturalCode.setGMS_Mr(stCalc.getGMS());
        //4: pKa
        constants.nextMode();
        naturalCode.setGMS_pKa(stCalc.getGMS());
        //5: pKb
        constants.nextMode();
        naturalCode.setGMS_pKb(stCalc.getGMS());
        // 6: pI
        constants.nextMode();
        naturalCode.setGMS_PI(stCalc.getGMS());
        this.naturalCode = naturalCode;
    }

    @Override
    Thread getThreadedCalculator() {
        return new Thread(() -> {
            GeneCode g=new GeneCode();
            Constants_Object constants = new Constants_Object();
            StabilityCalculator S=new StabilityCalculator(g, constants);
            while (true){
                int currentCode=getNextValue();
                if(currentCode>=codesWithoutValues.size())break;
                String[] rCode=codesWithoutValues.get(currentCode);
                MultiCharacteristicsCode codeData = new MultiCharacteristicsCode(rCode);
                g.changeCode(rCode);
                constants.reset();

                //0: polarReq
                codeData.setGMS_Polar(S.getGMS());
                //1: hydropathy
                constants.nextMode();
                codeData.setGMS_Hydro(S.getGMS());
                //2: molVol
                constants.nextMode();
                codeData.setGMS_MolVol(S.getGMS());
                //3: molWeight
                constants.nextMode();
                codeData.setGMS_Mr(S.getGMS());
                //4: pKa
                constants.nextMode();
                codeData.setGMS_pKa(S.getGMS());
                //5: pKb
                constants.nextMode();
                codeData.setGMS_pKb(S.getGMS());
                // 6: pI
                constants.nextMode();
                codeData.setGMS_PI(S.getGMS());

                codeData.calculateGreedyScore((MultiCharacteristicsCode)naturalCode);

                ValueBuffer[currentCode] = codeData;
            }
        });
    }

    void writeValueLine(FileWriter writer, int Number, IGreedyCode c){
        //Writes Value Line To File
        MultiCharacteristicsCode code = (MultiCharacteristicsCode) c;
        try {
            writer.write(Number+" "+String.join("~", c.getCode())+ " Score: "+c.getGreedy_score() +"["
                    + ToolMethods.df_short.format(code.getGMS_Polar()) +"|"
                    + ToolMethods.df_short.format(code.getGMS_Hydro()) +"|"
                    + ToolMethods.df_short.format(code.getGMS_MolVol()) +"|"
                    + ToolMethods.df_short.format(code.getGMS_Mr()) +"|"
                    + ToolMethods.df_short.format(code.getGMS_pKa()) +"|"
                    + ToolMethods.df_short.format(code.getGMS_pKb()) +"|"
                    + ToolMethods.df_short.format(code.getGMS_PI()) +"]\n");
        } catch (IOException e) {
            System.out.println("FileWriter Error");
        }
    }
}
