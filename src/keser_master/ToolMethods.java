package keser_master;

import Objects.Constants;
import org.biojava.nbio.core.sequence.DNASequence;

import java.text.DecimalFormat;
import java.util.List;

public class ToolMethods {
	public static DecimalFormat df = new DecimalFormat("0.0000");
	public static DecimalFormat df_long = new DecimalFormat("0.0000000");

	//calculates a value which determines how often a muatation results in a stop codon for single nucleotide mutations
	public static double[] getWeightedCountOfStopCodons_SNP(){
		double[] sum = new double[3];
		for(int a=0;a<4;a++){
			for(int b=0;b<4;b++){
				for(int c=0;c<4;c++){
					if(Constants.isStopCodon(a,b,c))continue;
					for(int x=0;x<4;x++){
						for (int pos=0;pos<3;pos++){
							String NewCodon= "";
							int basePrev;
							if (pos==0){
								NewCodon=Constants.Bases[x]+Constants.Bases[b]+Constants.Bases[c];
								basePrev=a;
							}else if (pos==1){
								NewCodon=Constants.Bases[a]+Constants.Bases[x]+Constants.Bases[c];
								basePrev=b;
							}else{
								NewCodon=Constants.Bases[a]+Constants.Bases[b]+Constants.Bases[x];
								basePrev=c;
							}
							if(!Constants.isStopCodon(NewCodon))continue;
							double count = 1;
							if (MainClass.nucleotideApriori!=null){
								count = count * MainClass.nucleotideApriori.getValue(basePrev);
							}
							if (MainClass.tripletApriori!=null){
								count = count * MainClass.tripletApriori.getValue(a,b,c);
							}
							sum[pos]=sum[pos]+count;
						}
					}
				}
			}
		}
		return sum;
	}

	//calculates a value which determines how often a muatation results in a stop codon for frameshift muatations
	public static double[] getWeightedCountOfStopCodons_Shift(){
		double[] sum= new double[2];
		for(int a=0;a<4;a++){
			for(int b=0;b<4;b++){
				for(int c=0;c<4;c++){
					if(Constants.isStopCodon(a,b,c))continue;
					for(int x=0;x<4;x++){
						for (int dir=0;dir<2;dir++){
							String NewCodon;
							if (dir==0){
								NewCodon=Constants.Bases[x]+Constants.Bases[a]+Constants.Bases[b];
							}else{
								NewCodon=Constants.Bases[b]+Constants.Bases[c]+Constants.Bases[x];
							}
							if(!Constants.isStopCodon(NewCodon))continue;
							double count = 1;
							if (MainClass.nucleotideTransition!=null){
								if (dir == 0){
									count = count * MainClass.nucleotideTransition.getValue(x,a);
								}else{
									count = count * MainClass.nucleotideTransition.getValue(c,x);
								}
							}
							if (MainClass.tripletTransition!=null){
								if (dir == 0){
									count = count * MainClass.tripletTransition.getValueLeft(a,b,c,x);
								}else{
									count = count * MainClass.tripletTransition.getValueRight(a,b,c,x);
								}

							}
							sum[dir]=sum[dir]+count;
						}
					}
				}
			}
		}
		return sum;
	}

	public static double getWeightedStopCodonFrequency_Overall(){
		double[] SNP=getWeightedCountOfStopCodons_SNP();
		double[] Shift= getWeightedCountOfStopCodons_Shift();
		return SNP[0]+SNP[1]+SNP[2]+Shift[0]+Shift[1];
	}

	//calculates the ratio of coding tripletts to stop-tripletts
	public static double getStopCodonFrequency(List<DNASequence> sequences, boolean inReadingFrame, int offset){
		int triplettcount = 0;
		int stopcodoncount = 0;
		for (DNASequence seq:sequences){
			String sequence = seq.getSequenceAsString();
			if(offset!=0){
				sequence=sequence.substring(offset);
			}
			int stepsize = 1;
			if(inReadingFrame){
				stepsize=3;
			}
			//count stopcodons
			for (int i = 0; i<sequence.length()-2; i = i + stepsize){
				String codon = sequence.substring(i,i+3);
				triplettcount++;
				if(Constants.isStopCodon(codon)){
					stopcodoncount++;
				}
			}

		}
		return (double)stopcodoncount/(double)triplettcount;
	}

	//Elementwise Difference between 2 4x4 matrices
	public static double[][] MatrixDiff(double[][] M1,double[][] M2){
		double[][]Erg=new double[4][4];
		for (int i=0;i<4;i++){
			for (int j=0;j<4;j++){
				Erg[i][j]=M1[i][j]-M2[i][j];
			}
		}
		double sum=0;
		for (int i=0;i<4;i++){
			for (int j=0;j<4;j++){
				sum+=Math.abs(Erg[i][j]);
			}
		}
		System.out.println("Gesamtdiferenz: "+df.format(sum));
		PrintMatrix(Erg);
		return Erg;
	}
	//Prints the 4x4 transition matrix in the console
	public static void PrintMatrix(double[][] M){
		System.out.println("Vertikal: s(n) horizontal: s(n+1)");
		System.out.println("--- T ------- C ------- A ------- G");
		System.out.println("T "+df.format(M[0][0])+" -- "+df.format(M[1][0])+" -- "+df.format(M[2][0])+" -- "+df.format(M[3][0]));
		System.out.println("C "+df.format(M[0][1])+" -- "+df.format(M[1][1])+" -- "+df.format(M[2][1])+" -- "+df.format(M[3][1]));
		System.out.println("A "+df.format(M[0][2])+" -- "+df.format(M[1][2])+" -- "+df.format(M[2][2])+" -- "+df.format(M[3][2]));
		System.out.println("G "+df.format(M[0][3])+" -- "+df.format(M[1][3])+" -- "+df.format(M[2][3])+" -- "+df.format(M[3][3]));	
	}
	//Prints a 4x4x4 matrix in the console
	public static void PrintMatrix(double[][][] M){
		String[] N={"T","C","A","G"};
		for (int a=0;a<4;a++){
			for (int b=0;b<4;b++){
				System.out.println(N[a]+N[b]+N[0]+": "+df.format(M[a][b][0])+" "+N[a]+N[b]+N[1]+": "+df.format(M[a][b][1])+" "+N[a]+N[b]+N[2]+": "+df.format(M[a][b][2])+" "+N[a]+N[b]+N[3]+": "+df.format(M[a][b][3])+" ");
			}
		}
	}
	//Prints NA Weightings
	public static void PrintMatrix(double[] M){
		System.out.println("T: "+M[0]);
		System.out.println("C: "+M[1]);
		System.out.println("A: "+M[2]);
		System.out.println("G: "+M[3]);
	}




}
