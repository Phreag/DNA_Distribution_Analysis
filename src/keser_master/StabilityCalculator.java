package keser_master;

import Objects.Constants;
import Objects.GeneCode;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;

public class StabilityCalculator {
	private GeneCode Code;
	private String[]Bases={"T","C","A","G"};
	
	private boolean baseAprioriWeighting=false;
	private double[] baseAprioriWeights;
	
	private boolean tripletaPrioriWeighting=false;
	private double[][][] tripletAprioriWeights;
	
	private boolean baseTransitionWeighting=false;
	private double[][] baseTransitionWeights;
	
	private boolean tripletTransitionWeighting=false;
	private double[][][][][]tripletTransitionWeights;

    private boolean printHistogram;

    public void setPrintHistogram(boolean printHistogram) {
        this.printHistogram = printHistogram;
    }

	//TransitionTransversionBias
	private int Bias=1;
	/* Deviationmode:
	 * 1=MS1 - Default
	 * 2=MS2
	 * 3=MS3
	 */
	public StabilityCalculator(GeneCode Code){
		this.Code=Code;
		this.Bias=MainClass.TransitionTransversionBias;
		//Changes the a Priori Weights for single Bases
		if (MainClass.baseAprioriEnabled){
			baseAprioriWeighting=true;
			baseAprioriWeights=MainClass.baseAprioriWeights;
		}
		//Changes the a Priori Weights for Triplets
		if (MainClass.tripletAprioriEnabled){
			tripletaPrioriWeighting=true;
			tripletAprioriWeights=MainClass.tripletAprioriWeights;
		}
		//Changes the Base-Transition Matrix
		if (MainClass.baseTransitionEnabled){
			baseTransitionWeighting=true;
			baseTransitionWeights=MainClass.baseTransitionWeights;
		}
		//Changes the Triplet-Transition Matrix used for Shift Calculations
		if (MainClass.tripletTransitionEnabled){
			tripletTransitionWeighting=true;
			tripletTransitionWeights=MainClass.tripletTransitionWeights;
		}
	}


	//Changes the Code used for Calculations
	public void ChangeCode(String[] Mapping){
		Code.changeCode(Mapping);
	}

	//Returns Deviations calculated with all current set Parameters
	//Only for single base Mutations
	//1=MS1
	//2=MS2
	//3=MS3
	public double get_BaseDeviation(int Modus){
		if(!(Modus>=1&&Modus<=3)){
			System.out.println("Bad Mode, Allowed: 1-3");
			return 0;
		}
		double deviation=0.0;
		for (int i=0;i<4;i++){
			String a=Bases[i];
			for (int j=0;j<4;j++){
				String b=Bases[j];
				for (int k=0;k<4;k++){
					String c=Bases[k];
					String codon1=a+b+c;
					String Amino=Code.getAminoAcid(codon1);
					if (Amino.length()!=3)continue; //Filters Stop Codons
					double Polar1=Constants.getPolarReq(Amino);
					Double Diff=0.0;
					for (int m=0;m<4;m++){
						String x=Bases[m];
						String Amino2="";
						String codon2="";
						switch (Modus){
						case 1:
						    codon2=x+b+c;
							break;
						case 2:
                            codon2=a+x+c;
							break;
						case 3:
						    codon2=a+b+x;
							break;
						}
                        Amino2=Code.getAminoAcid(codon2);
						double difference;
						if (Amino2.length()!=3){ //Filters Stop Codons or Applies NonsenseMutationWeighting

							continue;
							 
						}else{
							
							double Polar2=Constants.getPolarReq(Amino2);
							difference=(Polar1-Polar2)*(Polar1-Polar2);
						}
						
						if (baseAprioriWeighting){
							difference=difference*baseAprioriWeights[m];
						}
						if (tripletaPrioriWeighting){
							difference=difference*tripletAprioriWeights[i][j][k];
						}
						if(Bias!=1){
							String from="";
							switch (Modus){
							case 1:
								from=a;
								break;
							case 2:
								from=b;
								break;
							case 3:
								from=c;
								break;
							}
							if (isTransition(from, x)){
								difference=difference*Bias;
							}
						}
                        if(printHistogram){
                            printHiostogramEntry(difference,"SNP;"+codon1+"-"+codon2+";"+Amino+"-"+Amino2);
                        }
						Diff+=difference;
					}

					deviation=deviation+Diff;
				}
			}
		}

		switch (Modus){
		case 1:
			//61 codes * 3 = 183 - 9 Codes which can be mutated to a stop codon = 174
			//58 Mutations of these can be Transition
			deviation=deviation/(174+((Bias-1)*58));
			break;
		case 2:
			//61 codes * 3 = 183 - 7 Codes which can be mutated to a stop codon = 176
			//60 Mutations of these can be Transition
			deviation=deviation/(176+((Bias-1)*60));
			break;
		case 3:
			//61 codes * 3 = 183 - 7 Codes which can be mutated to a stop codon = 176
			//60 Mutations of these can be Transition
			deviation=deviation/(176+((Bias-1)*60));
			break;
		}
		return deviation;

		
	}
	//returns deviation for Shift Mutations
	//1=Right (+1), 2=Left (-1)
	public double get_ShiftDeviation(int Modus){
		if(!(Modus>=1&&Modus<=2)){
			System.out.println("Bad Mode, Allowed: 1-2");
			return 0;
		}
		double deviation=0.0;
		for (int i=0;i<4;i++){
			String a=Bases[i];
			for (int j=0;j<4;j++){
				String b=Bases[j];
				for (int k=0;k<4;k++){
					String c=Bases[k];
					String codon1=a+b+c;
					String Amino=Code.getAminoAcid(codon1);
					if (Amino.length()!=3)continue; //Filters Stop Codons as Origin
					double Polar1=Constants.getPolarReq(Amino);
					Double Diff=0.0;
					for (int m=0;m<4;m++){
						String x=Bases[m];
						String Amino2="";
						String codon2="";
						switch (Modus){
						case 1:
						    codon2=b+c+x;
							break;
						case 2:
                            codon2=x+a+b;
							break;
						}
                        Amino2=Code.getAminoAcid(codon2);
						double difference;
						if (Amino2.length()!=3){ //Filters Stop codons or applies NonsenseMutationWeights

							continue;
						}else{
							double Polar2=Constants.getPolarReq(Amino2);
							difference=(Polar1-Polar2)*(Polar1-Polar2);
						}
						if (baseAprioriWeighting){
							difference=difference*baseAprioriWeights[m];
						}
						if (tripletaPrioriWeighting){
							difference=difference*tripletAprioriWeights[i][j][k];
						}
						if (baseTransitionWeighting){
							switch(Modus){
							case 1:
								difference=difference*baseTransitionWeights[k][m];
								break;
							case 2:
								difference=difference*baseTransitionWeights[m][i];
								break;
							}
						}
						if(tripletTransitionWeighting){
							switch(Modus){
							case 1:
								difference=difference*tripletTransitionWeights[i][j][k][m][1];
								break;
							case 2:
								difference=difference*tripletTransitionWeights[i][j][k][m][0];
								break;
							}
						}
						if(printHistogram){
                            printHiostogramEntry(difference,"SHIFT;"+codon1+"-"+codon2+";"+Amino+"-"+Amino2);
                        }
						Diff+=difference;
					}
					deviation=deviation+Diff;
				}
			}
		}

		//244 Possible Mutations - 12 wich can result in a stop codon
		deviation=deviation/232;
		return deviation;
	}
	
	//returns true if the Mutation was a Transition
	private boolean isTransition(String from, String to){
		if(from.equalsIgnoreCase("A")&&to.equalsIgnoreCase("G"))return true;
		if(from.equalsIgnoreCase("G")&&to.equalsIgnoreCase("A"))return true;
		if(from.equalsIgnoreCase("C")&&to.equalsIgnoreCase("T"))return true;
		if(from.equalsIgnoreCase("T")&&to.equalsIgnoreCase("C"))return true;
		return false;
	}
	
	public double getWMS0(int Bias, double WMS1, double WMS2, double WMS3){
		WMS1=WMS1*(174+(Bias*58));
		WMS2=WMS2*(176+(Bias*60));
		WMS3=WMS3*(176+(Bias*60));
		return (WMS1+WMS2+WMS3)/((174+(Bias*58))+176+(Bias*60)+176+(Bias*60));
	}
	
	public double getMS0(double MS1, double MS2, double MS3){
		MS1=MS1*174;
		MS2=MS2*176;
		MS3=MS3*176;
		return(MS1+MS2+MS3)/(174+176+176);
	}
	
	public double getfMS(double rMS, double lMS){
		rMS=rMS*232;
		lMS=lMS*232;
		return(rMS+lMS)/(232+232);
	}
	public double getGMS(double MS1, double MS2, double MS3, double rMS, double lMS){
		MS1=MS1*174;
		MS2=MS2*176;
		MS3=MS3*176;
		rMS=rMS*232;
		lMS=lMS*232;
		return(MS1+MS2+MS3+rMS+lMS)/(174+176+176+232+232);
	}
	private void printHiostogramEntry(double value, String identifier){
        try {
            String Configuration="";
            if(MainClass.baseAprioriEnabled)Configuration=Configuration+"[NA]";
            if(MainClass.tripletAprioriEnabled)Configuration=Configuration+"[TA]";
            if(MainClass.baseTransitionEnabled)Configuration=Configuration+"[NT]";
            if(MainClass.tripletTransitionEnabled)Configuration=Configuration+"[TT]";
            if(Configuration.equals(""))Configuration="[NONE]";
            FileWriter fw = new FileWriter(new File("data/HIST_Detail_"+Configuration+".csv"), true);
            FileWriter fw2 = new FileWriter(new File("data/HIST_"+Configuration+".csv"), true);
            fw.write(String.valueOf(MainClass.df.format(value))+";"+identifier+"\n");
            fw2.write(String.valueOf(MainClass.df.format(value))+"\n");
            fw.close();
            fw2.close();
        } catch (IOException e) {
            System.out.println("Filewriter Error");
            e.printStackTrace();
        }
    }


}
