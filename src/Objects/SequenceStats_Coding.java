package Objects;

import java.text.DecimalFormat;

public class SequenceStats_Coding {
	private String Sequence;
	//Array Params: Base1. Base2, Base3, x, 0/1
	//x=Base in front or after
	//0= Base in Front of Triplet
	//1=Base after triplet
	private int[][][][][] rawData=new int[4][4][4][4][2];
	private String[]Bases={"T","C","A","G"};
	static DecimalFormat df = new DecimalFormat("0.0000"); 
	private double[] Base_aPriori;
	private double[][][] Triplet_aPriori;
	//NonsenseMutation Weights: Average Number of Tripletts following the current triplet
	private double[][][] NonsenseMutationWeights;
	private double [][][] NonsenseMutationWeights_raw =new double[4][4][4];
	private double[][][][][] TripletTransition;
	private double[][] BaseTransition;
	private int NumberOfTripletts;
	
	public SequenceStats_Coding (){
		
	}
	public void processSequence(String Sequence){
		this.Sequence=Sequence;
		NumberOfTripletts=Sequence.length()/3;
		updateData();
	}
	public void finalizeResults(){
		calculateBase_aPriori();
		calculateTriplet_aPriori();
		calculateBaseTransition();
		calculateTripletTransition();
		CalculateNonsenseMutationWeights();
		System.out.println("Finished Analyzing Coding Sequences");
	}

	private void updateData(){
		//Iterate over whole Sequence, Step Size 3;
		for (int i=0;i<Sequence.length()-3;i=i+3){
			try{
				int a=-1;
				int b=-1;
				int c=-1;
				int front=-1;
				int after=-1;
				
				// First Base of Triplet
				switch (Sequence.charAt(i)){
				case 'T':
					a=0;
					break;
				case 'C':
					a=1;
					break;
				case 'A':
					a=2;
					break;
				case 'G':
					a=3;
					break;
				default: 
					continue;
				}
				
				// Second Base of Triplet
				switch (Sequence.charAt(i+1)){
				case 'T':
					b=0;
					break;
				case 'C':
					b=1;
					break;
				case 'A':
					b=2;
					break;
				case 'G':
					b=3;
					break;
				default: 
					continue;
				}
				
				// Third Base of Triplet
				switch (Sequence.charAt(i+2)){
				case 'T':
					c=0;
					break;
				case 'C':
					c=1;
					break;
				case 'A':
					c=2;
					break;
				case 'G':
					c=3;
					break;
				default: 
					continue;
				}
				
				//Base After Triplet
				switch (Sequence.charAt(i+3)){
				case 'T':
					after=0;
					break;
				case 'C':
					after=1;
					break;
				case 'A':
					after=2;
					break;
				case 'G':
					after=3;
					break;
				}
				
				//Base Before Triplet
				if(i!=0){
					switch (Sequence.charAt(i-1)){
					case 'T':
						front=0;
						break;
					case 'C':
						front=1;
						break;
					case 'A':
						front=2;
						break;
					case 'G':
						front=3;
						break;
					}
				}
				int TriplettsAfterCurrentPos=NumberOfTripletts-(i/3+1);
				if(front==-1 && after==-1)continue;
				NonsenseMutationWeights_raw[a][b][c]+=TriplettsAfterCurrentPos;
				if (front!=-1){
					rawData[a][b][c][front][0]++;
				}
				if (after!=-1){
					rawData[a][b][c][after][1]++;
				}
			}catch(Exception e){}
		}
	}
	
	private void CalculateNonsenseMutationWeights(){
		//Average on each Triplet
		//Normalize to 1 (use as Factor)
		int totalsum=0;
		int totalcount=0;
		NonsenseMutationWeights=new double[4][4][4];
		for(int i=0;i<4;i++){
			for(int j=0;j<4;j++){
				for(int k=0;k<4;k++){
					int count=rawData[i][j][k][0][1]+rawData[i][j][k][1][1]+rawData[i][j][k][2][1]+rawData[i][j][k][3][1];
					totalcount+=count;
					NonsenseMutationWeights[i][j][k]=NonsenseMutationWeights_raw[i][j][k]/count;
					totalsum+=NonsenseMutationWeights_raw[i][j][k];
				}
			}
		}
		double normalizedsum=0;
		for(int i=0;i<4;i++){
			for(int j=0;j<4;j++){
				for(int k=0;k<4;k++){
					normalizedsum= normalizedsum+NonsenseMutationWeights[i][j][k];
				}
			}
		}
		double average = normalizedsum/64;
		System.out.println(average);
		for(int i=0;i<4;i++){
			for(int j=0;j<4;j++){
				for(int k=0;k<4;k++){
					NonsenseMutationWeights[i][j][k]=NonsenseMutationWeights[i][j][k]/average;
				}
			}
		}
	}
	
	//Calculates NA weightings - Edited for 3-Step
		private void calculateBase_aPriori(){
			int[] BaseCount=new int[4];
			int overallCount=0;
			Base_aPriori=new double[4];
			for(int i=0;i<4;i++){
				for(int j=0;j<4;j++){
					for(int k=0;k<4;k++){
						for(int l=0;l<4;l++){
							BaseCount[i]+=rawData[i][j][k][l][1];
							BaseCount[j]+=rawData[i][j][k][l][1];
							BaseCount[k]+=rawData[i][j][k][l][1];
							overallCount+=rawData[i][j][k][l][1]*3;
						}
					}
				}
			}
			for (int i=0;i<4;i++){
				Base_aPriori[i]=((double)BaseCount[i]/(double)overallCount)*4;
			}
			System.out.println("Base_aPriori Average: "+(Base_aPriori[0]+Base_aPriori[1]+Base_aPriori[2]+Base_aPriori[3])/4);
		}
		//Calculates TA weightings	
		private void calculateTriplet_aPriori(){
			int[][][] TripletCount=new int[4][4][4];
			int overallCount=0;
			Triplet_aPriori=new double[4][4][4];
			for(int i=0;i<4;i++){
				for(int j=0;j<4;j++){
					for(int k=0;k<4;k++){
						for(int l=0;l<4;l++){
							//This will ignore 4 Bases at the beginning - shit happens :P
							String Pattern=Bases[i]+Bases[j]+Bases[k];
							if (Pattern.equalsIgnoreCase("TAA")||Pattern.equalsIgnoreCase("TAG")||Pattern.equalsIgnoreCase("TGA"))continue;
							TripletCount[i][j][k]+=rawData[i][j][k][l][1];
							overallCount+=rawData[i][j][k][l][1];
						}
					}
				}
			}
			double sum=0;
			for(int i=0;i<4;i++){
				for(int j=0;j<4;j++){
					for(int k=0;k<4;k++){
						Triplet_aPriori[i][j][k]=((double)TripletCount[i][j][k]*61)/(double)overallCount;
						sum+=Triplet_aPriori[i][j][k];
					}
				}
			}
			sum=sum/61;
			System.out.println("Triplet_aPriori Average: "+sum);
		}
		//Calculates NT weightings - edited for 3-Step	
		private void calculateBaseTransition(){
			int[][] BaseTransitionCount=new int[4][4];
			int overallCount=0;
			BaseTransition=new double[4][4];
			for(int i=0;i<4;i++){
				for(int j=0;j<4;j++){
					for(int k=0;k<4;k++){
						for(int l=0;l<4;l++){
							BaseTransitionCount[i][j]+=rawData[i][j][k][l][1];
							BaseTransitionCount[j][k]+=rawData[i][j][k][l][1];
							BaseTransitionCount[k][l]+=rawData[i][j][k][l][1];
							overallCount+=rawData[i][j][k][l][1]*3;
						}
					}
				}
			}
			double sum=0;
			for (int i=0;i<4;i++){
				for (int j=0;j<4;j++){
					BaseTransition[i][j]=((double)BaseTransitionCount[i][j]/(double)overallCount)*16;
					sum+=BaseTransition[i][j];
				}
			}
			sum=sum/16;
			System.out.println("BaseTransition Average: "+sum);
		}
		//Calculates TT weightings	
		private void calculateTripletTransition(){
			int overallCountFront=0;
			int overallCountAfter=0;
			TripletTransition=new double[4][4][4][4][2];
			for(int i=0;i<4;i++){
				for(int j=0;j<4;j++){
					for(int k=0;k<4;k++){
						for(int l=0;l<4;l++){
							overallCountFront+=rawData[i][j][k][l][0];
							overallCountAfter+=rawData[i][j][k][l][1];
						}
					}
				}
			}
			double sumFront=0;
			double sumAfter=0;
			for(int i=0;i<4;i++){
				for(int j=0;j<4;j++){
					for(int k=0;k<4;k++){
						for(int l=0;l<4;l++){
							TripletTransition[i][j][k][l][0]=((double)rawData[i][j][k][l][0]/(double)overallCountFront)*256;
							TripletTransition[i][j][k][l][1]=((double)rawData[i][j][k][l][1]/(double)overallCountAfter)*256;
							sumFront+=TripletTransition[i][j][k][l][0];
							sumAfter+=TripletTransition[i][j][k][l][1];
						}
					}
				}
			}
			sumFront=sumFront/256;
			sumAfter=sumAfter/256;
			System.out.println("Triplet Transition Average: Front:"+sumFront+" After:"+sumAfter);
		}
	
	//Difference by Element between matrices
	

	public double[] getBase_aPriori() {
		return Base_aPriori;
	}

	public double[][][] getTriplet_aPriori() {
		return Triplet_aPriori;
	}

	public double[][][][][] getTripletTransition() {
		return TripletTransition;
	}

	public double[][] getBaseTransition() {
		return BaseTransition;
	}
	public double[][][] getNonsenseMutationWeights() {
		return NonsenseMutationWeights;
	}
}
