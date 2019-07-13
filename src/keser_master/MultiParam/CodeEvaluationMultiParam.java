package keser_master.MultiParam;

import keser_master.MainClass;
import keser_master.MultiParam.code.MultiCharacteristicsCode;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;

public class CodeEvaluationMultiParam {
    private MultiCharacteristicsCode[] values;
    boolean firstrun = true;
    int[] betterCodes = new int[8];
    MultiCharacteristicsCode naturalCode;
    int codesSeen = 0;

    public CodeEvaluationMultiParam(MultiCharacteristicsCode[] values) {
        this.values = values;
    }

    public void setValues(MultiCharacteristicsCode[] values) {
        this.values = values;
    }

    public int[] countBetterCodes(String Title) {
        System.out.println("Evaluating Results and counting better codes found...");

        //Natural code needs to be #1

        if (firstrun) {
            naturalCode = values[0];
            try {
                FileWriter fw = new FileWriter(new File("data/EvaluationResultsMultiParam.log"), true);
                fw.write("\n");
                fw.write("Calculation Start: " + new SimpleDateFormat("dd.MM.yyyy HH:mm:ss").format(new Date()) + "("+Title+")\n");
                fw.write("   Weightings:" + MainClass.getConfigString() + "\n");
                fw.close();
            } catch (IOException e) {
                System.out.println("Filewriter Error");
                e.printStackTrace();
            }
            printCodeData(naturalCode, false);
        }

        for (int i = firstrun ? 0 : 1; i < values.length; i++) {
            boolean isCompleteBetter = true;

            //0: polarReq
            if (values[i].getGMS_Polar() < naturalCode.getGMS_Polar()) {
                betterCodes[0]++;
            } else {
                isCompleteBetter = false;
            }
            //1: hydropathy
            if (values[i].getGMS_Hydro() < naturalCode.getGMS_Hydro()) {
                betterCodes[1]++;
            } else {
                isCompleteBetter = false;
            }
            //2: molVol
            if (values[i].getGMS_MolVol() < naturalCode.getGMS_MolVol()) {
                betterCodes[2]++;
            } else {
                isCompleteBetter = false;
            }
            //3: molWeight
            if (values[i].getGMS_Mr() < naturalCode.getGMS_Mr()) {
                betterCodes[3]++;
            } else {
                isCompleteBetter = false;
            }
            //4: pKa
            if (values[i].getGMS_pKa() < naturalCode.getGMS_pKa()) {
                betterCodes[4]++;
            } else {
                isCompleteBetter = false;
            }
            //5: pKb
            if (values[i].getGMS_pKb() < naturalCode.getGMS_pKb()) {
                betterCodes[5]++;
            } else {
                isCompleteBetter = false;
            }
            // 6: pI
            if (values[i].getGMS_PI() < naturalCode.getGMS_PI()) {
                betterCodes[6]++;
            } else {
                isCompleteBetter = false;
            }


            if (isCompleteBetter) {
                betterCodes[7]++;
                printCodeData(values[i], true);
            }
            if(firstrun || i != 0){
                codesSeen++;
            }

            if(codesSeen % 100000 == 0){
                System.out.println(new Date().toString()+": Codes Analyzed: "+codesSeen+" BetterCodes: "+Arrays.toString(betterCodes));
            }
        }

        firstrun = false;
        return betterCodes;
    }

    public void printSummary(){
        try {
            FileWriter fw = new FileWriter(new File("data/EvaluationResultsMultiParam.log"), true);
            fw.write("######## SUMMARY ########\n");
            fw.write("   Date | #CodesSeen | [polarReq, hydropathy, molVol, molWeight, pKa, pKb, pI, ImmerBesser]\n");
            fw.write("   " + new SimpleDateFormat("dd.MM.yyyy HH:mm:ss").format(new Date())+ " | " +codesSeen+ " | " + Arrays.toString(betterCodes) + "\n");
            fw.close();
        } catch (IOException e) {
            System.out.println("Filewriter Error");
            e.printStackTrace();
        }
    }

    private void printCodeData(MultiCharacteristicsCode code, boolean betterCode) {
        try {
            FileWriter fw = new FileWriter(new File("data/EvaluationResultsMultiParam.log"), true);
            if (betterCode) {
                fw.write("!!! Better Code: !!!\n");
            } else {
                fw.write("   Natural Code:\n");
            }
            fw.write("   Mapping: " + Arrays.toString(code.getCode()) + "\n");
            fw.write("   GMS Scores [polarReq, hydropathy, molVol, molWeight, pKa, pKb, pI]\n");
            fw.write("   [" + code.getGMS_Polar() +
                    ", " + code.getGMS_Hydro() +
                    ", " + code.getGMS_MolVol() +
                    ", " + code.getGMS_Mr() +
                    ", " + code.getGMS_pKa() +
                    ", " + code.getGMS_pKb() +
                    ", " + code.getGMS_PI() + "]\n");

            fw.write("\n");
            fw.close();
        } catch (IOException e) {
            System.out.println("Filewriter Error");
            e.printStackTrace();
        }
    }


}
