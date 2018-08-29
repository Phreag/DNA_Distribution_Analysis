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

import Objects.*;

public class MainClass {
	static String[] Code={"Leu","Pro","His","Gln","Arg","Ile","Met","Thr","Asn","Lys","Ser","Val","Ala","Asp","Glu","Gly","Phe","Tyr","Cys","Trp"};
	static DecimalFormat df = new DecimalFormat("0.0000"); 
	static GenBankConnection conn=new GenBankConnection();
	public static double[] baseAprioriWeights;
	public static double[][][] tripletAprioriWeights;
	public static double[][][] nonsenseMutationWeights;
	public static double[][] baseTransitionWeights;
	public static double[][][][][]tripletTransitionWeights;
	public static boolean baseAprioriEnabled=false;
	public static boolean tripletAprioriEnabled=false;
	public static boolean baseTransitionEnabled=false;
	public static boolean tripletTransitionEnabled=false;
	public static boolean nonsenseWeightingEnabled=false;
	public static int TransitionTransversionBias=1;
	//Nonsense Mutation Factor: Sets weight for the Error produced by Nonsense Mutation
	//as Factor * NumberOfTripletsAfter Sequence
	public static double NonsenseMutationFactor=1;
	//This File Should be run with at least 7Gb of Java Heap Space!
	public static void main (String[] args){
		//###################################################################
		//#########################   DEBUG AREA   ##########################
		//###################################################################
		/*
		 * 
		 * ToDo:
		 * Nucleotidverteilungen auf Gesamter Sequenz vs Codierende Sequenz: Tabelle 3.2
		 */
		//CCDS
//		List<DNASequence> cDNA=conn.LoadMixedFile();
//		System.out.println("Size: " +cDNA.size());
//		SequenceStats_Coding Stat=new SequenceStats_Coding();
//		for (DNASequence Seq : cDNA){
//			Stat.ProcessSequence(Seq.getSequenceAsString());
//		}
//		Stat.FinalizeResults();
//		System.out.println("CCDS NA Gewichte");
//		Stat.PrintMatrix(Stat.getBase_aPriori());
//		//Chromosom 1
//		DNASequence Seq1=conn.LoadFastaFile("NC_000001.11");
//		SequenceStats Stat2=new SequenceStats(Seq1.getSequenceAsString());
//		System.out.println("Chromosom 1 NA Gewichte");
//		Stat.PrintMatrix(Stat2.getBase_aPriori());
		/*Resultat
		 * CCDS NA Gewichte T: 0.8669436165203899 C: 1.0655803531259607 A: 1.0047852038259497 G: 1.0626908265276998
		 * Chromosom 1 NA Gewichte T: 1.1670230633493237 C: 0.8339956930420922 A: 1.1640052355203168 G: 0.8349760080882673
		 */
		
		 /*
		 * Entstehen Statistisch gesehen durch Gewichtungen mehr Stoppcondons?
		 * Hypothese: Die DNA ist Optimiert dass es schnell zum Abbruch kommt
		 * Vergleich: Komplette Sequenz vs Codierende Regionen
		 */
		//Stoppcodons je Base bei Punktmutation Chromosom 1
//		DNASequence Seq1=conn.LoadFastaFile("NC_000001.11");
//		SequenceStats Stat=new SequenceStats(Seq1.getSequenceAsString());
//		System.out.println("Chromosom 1: Keine Gewichtungen");
//		ToolMethods.getWeightedCountOfStopCodons();
//		baseAprioriWeights=Stat.getBase_aPriori();
//		tripletAprioriWeights=Stat.getTriplet_aPriori();
//		setWeightings(true,false,false,false,false);
//		System.out.println("Chromosom 1: NA");
//		ToolMethods.getWeightedCountOfStopCodons();
//		setWeightings(false,true,false,false,false);
//		System.out.println("Chromosom 1: TA");
//		ToolMethods.getWeightedCountOfStopCodons();
//		setWeightings(true,true,false,false,false);
//		System.out.println("Chromosom 1: NA+TA");
//		ToolMethods.getWeightedCountOfStopCodons();
		/*
		 * Chromosom 1: Keine Gewichtungen
		 * Pos1: 9,0000 Pos2: 7,0000 Pos3: 7,0000
		 * Chromosom 1: NA
		 * Pos1: 8,4989 Pos2: 6,8380 Pos3: 6,8380
		 * Chromosom 1: TA
		 * Pos1: 11,2169 Pos2: 7,3554 Pos3: 7,6613
		 * Chromosom 1: NA+TA
		 * Pos1: 11,0398 Pos2: 7,3938 Pos3: 7,6449

		 */
		//Stoppcodons je Base bei Punktmutation CCDS
		List<DNASequence> cDNA=conn.LoadMixedFile();
		System.out.println("Size: " +cDNA.size());
		SequenceStats_Coding Stat=new SequenceStats_Coding();
		for (DNASequence Seq : cDNA){
			Stat.ProcessSequence(Seq.getSequenceAsString());
		}
		Stat.FinalizeResults();
		System.out.println("CCDS: Keine Gewichtungen");
		ToolMethods.getWeightedCountOfStopCodons();
		baseAprioriWeights=Stat.getBase_aPriori();
		tripletAprioriWeights=Stat.getTriplet_aPriori();
		setWeightings(true,false,false,false,false);
		System.out.println("CCDS: NA");
		ToolMethods.getWeightedCountOfStopCodons();
		setWeightings(false,true,false,false,false);
		System.out.println("CCDS: TA");
		ToolMethods.getWeightedCountOfStopCodons();
		setWeightings(true,true,false,false,false);
		System.out.println("CCDS: NA+TA");
		ToolMethods.getWeightedCountOfStopCodons();
		 /* Vergleich: TT Gewichtung Nat. Code vs Random Set
		 * Bei Shiftmutation. Keine Gewichtung: 267 besser
		 * TT komplette Sequenz: 28 besser
		 * TT CCDS: lms 1502 rms 1191 fms 965
		 */
		
		 /* Tabelle 3.11 reproduzieren
		 */
		
		 /* WMS0 Vergleich abhängig von Transition/Transversion Bias, Tabelle 3.12
		 */
		//Analysis of Coding Sequences in Reading frame
//		List<DNASequence> cDNA=conn.LoadMixedFile();
//		System.out.println("Size: " +cDNA.size());
//		SequenceStats_Coding Stat=new SequenceStats_Coding();
//		for (DNASequence Seq : cDNA){
//			Stat.ProcessSequence(Seq.getSequenceAsString());
//		}
//		Stat.FinalizeResults();
//		
//		CodePermutation P=new CodePermutation();
//		P.loadDefaultcodeSet();
//		new CodeEvaluation(P.calculateValues()).countBetterCodes();
//		baseAprioriWeights=Stat.getBase_aPriori();
//		tripletAprioriWeights=Stat.getTriplet_aPriori();
//		baseTransitionWeights=Stat.getBaseTransition();
//    	tripletTransitionWeights=Stat.getTripletTransition();
// 		nonsenseMutationWeights=Stat.getNonsenseMutationWeights();
// 		NonsenseMutationFactor=100.0;
// 		setWeightings(false,false,false,false,true);
// 		P=new CodePermutation();
//		P.loadDefaultcodeSet();
//		new CodeEvaluation(P.calculateValues()).countBetterCodes();
		
//		CodePermutation P=new CodePermutation();
//		P.loadDefaultcodeSet();
//		//new CodeEvaluation(P.calculateValues()).countBetterCodes();
//		baseAprioriWeights=Stat.getBase_aPriori();
//		tripletAprioriWeights=Stat.getTriplet_aPriori();
//		baseTransitionWeights=Stat.getBaseTransition();
//		tripletTransitionWeights=Stat.getTripletTransition();
//		nonsenseMutationWeights=Stat.getNonsenseMutationWeights();
//		NonsenseMutationFactor=10.0;
//		setWeightings(false,false,false,false,true);			
//		new CodeEvaluation(P.calculateValues()).countBetterCodes();
				
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
//  		DNASequence Seq1=conn.LoadFastaFile("NC_000001.11");
//  		SequenceStats Stat=new SequenceStats(Seq1.getSequenceAsString());
//		
//		baseAprioriWeights=Stat.getBase_aPriori();
//		tripletAprioriWeights=Stat.getTriplet_aPriori();
//		baseTransitionWeights=Stat.getBaseTransition();
//		tripletTransitionWeights=Stat.getTripletTransition();
//		
//		String[]Bases={"T","C","A","G"};
//		DecimalFormat df = new DecimalFormat("0.0000"); 
//		for (int x=0;x<4;x++){
//			for (int y=0;y<4;y++){
//				for (int z=0;z<4;z++){
//						System.out.println(Bases[x]+Bases[y]+Bases[z]+";"+df.format(tripletTransitionWeights[x][y][z][0][0]/4)+";"+df.format(tripletTransitionWeights[x][y][z][0][1]/4)+";"+df.format(tripletTransitionWeights[x][y][z][1][0]/4)+";"+df.format(tripletTransitionWeights[x][y][z][1][1]/4)+";"+df.format(tripletTransitionWeights[x][y][z][2][0]/4)+";"+df.format(tripletTransitionWeights[x][y][z][2][1]/4)+";"+df.format(tripletTransitionWeights[x][y][z][3][0]/4)+";"+df.format(tripletTransitionWeights[x][y][z][3][1]/4));
//				}
//			}
//		}
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
	private static void setWeightings(boolean ba, boolean ta, boolean bt, boolean tt, boolean nm){
		baseAprioriEnabled=ba;
		tripletAprioriEnabled=ta;
		baseTransitionEnabled=bt;
		tripletTransitionEnabled=tt;
		nonsenseWeightingEnabled=nm;
	}
}
