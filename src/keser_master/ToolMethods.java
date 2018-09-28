package keser_master;

import java.text.DecimalFormat;

public class ToolMethods {
	static DecimalFormat df = new DecimalFormat("0.0000"); 
	private static String[]Bases={"T","C","A","G"};
	
	public static double[] getWeightedCountOfStopCodons_SNP(){
		double[] sum= new double[3];
		for(int a=0;a<4;a++){
			for(int b=0;b<4;b++){
				for(int c=0;c<4;c++){
					
					String Pattern=Bases[a]+Bases[b]+Bases[c];
					if (Pattern.equalsIgnoreCase("TAA")||Pattern.equalsIgnoreCase("TAG")||Pattern.equalsIgnoreCase("TGA"))continue;
					for(int x=0;x<4;x++){
						for (int pos=0;pos<3;pos++){
							String NewCodon= "";
							int basePrev;
							if (pos==0){
								NewCodon=Bases[x]+Bases[b]+Bases[c];
								basePrev=a;
							}else if (pos==1){
								NewCodon=Bases[a]+Bases[x]+Bases[c];
								basePrev=b;
							}else{
								NewCodon=Bases[a]+Bases[b]+Bases[x];
								basePrev=c;
							}
							//No Stopcodon
							if (!NewCodon.equalsIgnoreCase("TAA") && !NewCodon.equalsIgnoreCase("TAG") && !NewCodon.equalsIgnoreCase("TGA"))continue;
							double count = 1;
							if (MainClass.baseAprioriEnabled){
								count = count * MainClass.baseAprioriWeights[basePrev];
							}
							if (MainClass.tripletAprioriEnabled){
								count = count * MainClass.tripletAprioriWeights[a][b][c];
							}
							sum[pos]=sum[pos]+count;
						}
					}
				}
			}
		}
		return sum;
	}
	public static double[] getWeightedCountOfStopCodons_Shift(){
		double[] sum= new double[2];
		for(int a=0;a<4;a++){
			for(int b=0;b<4;b++){
				for(int c=0;c<4;c++){

					String Pattern=Bases[a]+Bases[b]+Bases[c];
					if (Pattern.equalsIgnoreCase("TAA")||Pattern.equalsIgnoreCase("TAG")||Pattern.equalsIgnoreCase("TGA"))continue;
					for(int x=0;x<4;x++){


						for (int dir=0;dir<2;dir++){
							String NewCodon= "";
							int basePrev;
							if (dir==0){
								NewCodon=Bases[x]+Bases[a]+Bases[b];
							}else{
								NewCodon=Bases[b]+Bases[c]+Bases[x];
							}
							//No Stopcodon
							if (!NewCodon.equalsIgnoreCase("TAA") && !NewCodon.equalsIgnoreCase("TAG") && !NewCodon.equalsIgnoreCase("TGA"))continue;
							double count = 1;
							if (MainClass.baseTransitionEnabled){
								if (dir == 0){
									count = count * MainClass.baseTransitionWeights[x][a];
								}else{
									count = count * MainClass.baseTransitionWeights[c][x];
								}
							}
							if (MainClass.tripletTransitionEnabled){
								if (dir == 0){
									count = count * MainClass.tripletTransitionWeights[a][b][c][x][0];
								}else{
									count = count * MainClass.tripletTransitionWeights[a][b][c][x][1];
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
		double total=SNP[0]+SNP[1]+SNP[2]+Shift[0]+Shift[1];
		return total;
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
