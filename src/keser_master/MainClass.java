package keser_master;
 
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.MalformedURLException;
import java.net.URL;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map.Entry;
import org.apache.commons.io.FileUtils;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.io.FastaReader;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.io.GenericFastaHeaderParser;
import org.biojava.nbio.core.sequence.io.ProteinSequenceCreator;

import Objects.GeneCode;
import Objects.SequenceStats;

public class MainClass {
	static String[] Code={"Leu","Pro","His","Gln","Arg","Ile","Met","Thr","Asn","Lys","Ser","Val","Ala","Asp","Glu","Gly","Phe","Tyr","Cys","Trp"};
	static DecimalFormat df = new DecimalFormat("0.0000"); 
	static GenBankConnection conn=new GenBankConnection();
	public static double[] baseAprioriWeights;
	public static double[][][] tripletAprioriWeights;
	public static double[][] baseTransitionWeights;
	public static double[][][][][]tripletTransitionWeights;
	public static boolean baseAprioriEnabled=false;
	public static boolean tripletAprioriEnabled=false;
	public static boolean baseTransitionEnabled=false;
	public static boolean tripletTransitionEnabled=false;
	public static int TransitionTransversionBias=1;
	//This File Should be run with at least 7Gb of Java Heap Space!
	public static void main (String[] args){
		//###################################################################
		//#########################   DEBUG AREA   ##########################
		//###################################################################
		
		//Get One miexd String from Multi-Sequence File
		//		List<DNASequence> Mixed=conn.LoadMixedFile();
		//		StringBuilder builder = new StringBuilder();
		//		for(DNASequence seq : Mixed) {		    
		//		builder.append(seq.getSequenceAsString());
		//    	}
				
		//		String MixedSeq=builder.toString();
		//		System.out.println("Mixed Length: "+MixedSeq.length());		
		//		double[] w=stat.getNucleotideDistribution(MixedSeq);
		//		double[] tempfactors={w[0]/0.25,w[1]/0.25,w[2]/0.25,w[3]/0.25};
		//		factors=tempfactors;
		//		tweights=stat.getTripletDistribution(MixedSeq);
		
		//Batch Calculation to compare different sequences and their impact on the Scores
		//		List<String> Sequences=new ArrayList<String>();
		
		//		Sequences.add("NC_001806"); //Herpes
		//		Sequences.add("NC_001591"); //Papillomavirus
		//		Sequences.add("NC_003461"); //Parainfluenza
		//		Sequences.add("NC_001959"); //Norovirus G1
		//		Sequences.add("NC_002031");//Gelbfieber
		//		Sequences.add("NC_012532");//Zika
		//		Sequences.add("NC_001477");//Dengue
		//		Sequences.add("NC_004102");//Hepatitis C
		//		Sequences.add("NC_002549");//Zaire Ebola
		
				
		//		for (String Str:Sequences){
		//			DNASequence Seq1=conn.LoadFastaFile(Str);
		//			SequenceStats Stat=new SequenceStats(Seq1.getSequenceAsString());
		//			baseAprioriWeights=Stat.getBase_aPriori();
		//			tripletAprioriWeights=Stat.getTriplet_aPriori();
		//			baseTransitionWeights=Stat.getBaseTransition();
		//			tripletTransitionWeights=Stat.getTripletTransition();
		//			TransitionTransversionBias=2;
		//			setWeightings(true, true, true, true);
		//			GeneCode G=new GeneCode(Code);
		//			StabilityCalculator S=new StabilityCalculator(G);
		//			double MS1=S.get_BaseDeviation(1);//MS1
		//			double MS2=S.get_BaseDeviation(2);//MS2
		//			double MS3=S.get_BaseDeviation(3);//MS3
		//			double MS0=S.getMS0(MS1, MS2, MS3);//MS0
		//			double rMS=S.get_ShiftDeviation(1);//rMS
		//			double lMS=S.get_ShiftDeviation(2);//lMS
		//			double fMS=S.getfMS(rMS, lMS);//fMS
		//			double GMS=S.getGMS(MS1, MS2, MS3, rMS, lMS);
		//			FileWriter Log;
		//			try {
		//				Log = new FileWriter("data/Log.txt", true);
		//				Log.write(Str+", "+MS1+", "+MS2+", "+MS3+", "+MS0+", "+rMS+", "+lMS+", "+fMS+", "+GMS+"\n");
		//				Log.close();
		//			} catch (IOException e) {
		//				e.printStackTrace();
		//			}
		//		}
  		DNASequence Seq1=conn.LoadFastaFile("NC_000001.11");
  		SequenceStats Stat=new SequenceStats(Seq1.getSequenceAsString());
		
		baseAprioriWeights=Stat.getBase_aPriori();
		tripletAprioriWeights=Stat.getTriplet_aPriori();
		baseTransitionWeights=Stat.getBaseTransition();
		tripletTransitionWeights=Stat.getTripletTransition();
		
		String[]Bases={"T","C","A","G"};
		DecimalFormat df = new DecimalFormat("0.0000"); 
		for (int x=0;x<4;x++){
			for (int y=0;y<4;y++){
				for (int z=0;z<4;z++){
						System.out.println(Bases[x]+Bases[y]+Bases[z]+";"+df.format(tripletTransitionWeights[x][y][z][0][0]/4)+";"+df.format(tripletTransitionWeights[x][y][z][0][1]/4)+";"+df.format(tripletTransitionWeights[x][y][z][1][0]/4)+";"+df.format(tripletTransitionWeights[x][y][z][1][1]/4)+";"+df.format(tripletTransitionWeights[x][y][z][2][0]/4)+";"+df.format(tripletTransitionWeights[x][y][z][2][1]/4)+";"+df.format(tripletTransitionWeights[x][y][z][3][0]/4)+";"+df.format(tripletTransitionWeights[x][y][z][3][1]/4));
				}
			}
		}
//		TransitionTransversionBias=2;
////		CodePermutation P=new CodePermutation();
////		P.generateCodes();
//		setWeightings(false, true, false, true);
//		CodeFinder C=new CodeFinder();
//    	C.RunCodeFinder(15);
//  	new CodeEvaluation(P.calculateValues()).countBetterCodes();
//		CodeFinder C=new CodeFinder();
//      setWeightings(true,true,true,true);
//		CodePermutation P=new CodePermutation();
//		new CodeEvaluation(P.calculateValues()).countBetterCodes();
//		C.RunCodeFinder(15);
//		
//		CodePermutation P=new CodePermutation();
//		P.generateCodes();
//		
//  		CodeFinder C=new CodeFinder();
//  		C.RunCodeFinder(20);
//		//Ohne Gewichtung
//		CodePermutation P=new CodePermutation();
//		new CodeEvaluation(P.calculateValues()).countBetterCodes();
//		

		
		

		//###################################################################
		//###################################################################
		//###################################################################
		

	}
	//Set Weightings here to enable or disable them globally
	private static void setWeightings(boolean ba, boolean ta, boolean bt, boolean tt){
		baseAprioriEnabled=ba;
		tripletAprioriEnabled=ta;
		baseTransitionEnabled=bt;
		tripletTransitionEnabled=tt;
	}
	
}
