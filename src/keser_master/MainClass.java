package keser_master;
 
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.MalformedURLException;
import java.net.URL;
import java.sql.SQLOutput;
import java.text.DecimalFormat;
import java.util.Arrays;
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
	public static double[][] baseTransitionWeights;
	public static double[][][][][]tripletTransitionWeights;
	public static boolean baseAprioriEnabled=false;
	public static boolean tripletAprioriEnabled=false;
	public static boolean baseTransitionEnabled=false;
	public static boolean tripletTransitionEnabled=false;
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
		 * ToDo: Nucleotidverteilungen auf Gesamter Sequenz vs Codierende Sequenz: Tabelle 3.2
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
		 * ToDo: Entstehen Statistisch gesehen durch Gewichtungen mehr Stoppcondons?
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
//		List<DNASequence> cDNA=conn.LoadMixedFileReadingframe();
//		System.out.println("Size: " +cDNA.size());
//		SequenceStats_Coding Stat=new SequenceStats_Coding();
//		for (DNASequence Seq : cDNA){
//			Stat.ProcessSequence(Seq.getSequenceAsString());
//		}
//		Stat.FinalizeResults();
//		System.out.println("CCDS: Keine Gewichtungen");
//		ToolMethods.getWeightedCountOfStopCodons();
//		baseAprioriWeights=Stat.getBase_aPriori();
//		tripletAprioriWeights=Stat.getTriplet_aPriori();
//		setWeightings(true,false,false,false,false);
//		System.out.println("CCDS: NA");
//		ToolMethods.getWeightedCountOfStopCodons();
//		setWeightings(false,true,false,false,false);
//		System.out.println("CCDS: TA");
//		ToolMethods.getWeightedCountOfStopCodons();
//		setWeightings(true,true,false,false,false);
//		System.out.println("CCDS: NA+TA");
//		ToolMethods.getWeightedCountOfStopCodons();
		 /* Vergleich: TT Gewichtung Nat. Code vs Random Set
		 * Bei Shiftmutation. Keine Gewichtung: 267 besser
		 * TT komplette Sequenz: 28 besser
		 * TT CCDS: lms 1502 rms 1191 fms 965
		 */
		
		 /*ToDo:  Tabelle 3.11 reproduzieren. Je 1 mal mit Homo Sapiens, E.Coli und Ciona intestinalis
		 */
/*		List<DNASequence> cDNA=conn.LoadMixedFileReadingframe("Escherichia_coli.HUSEC2011CHR1.cdna.all.fasta");
		System.out.println("Size: " +cDNA.size());
		SequenceStats_Coding Stat=new SequenceStats_Coding();
		for (DNASequence Seq : cDNA){
			Stat.ProcessSequence(Seq.getSequenceAsString());
		}
		Stat.FinalizeResults();
        baseAprioriWeights=Stat.getBase_aPriori();
        tripletAprioriWeights=Stat.getTriplet_aPriori();
        tripletTransitionWeights=Stat.getTripletTransition();
		String T="Tabelle 3.11 Reproduktion: Escherichia_coli.HUSEC2011CHR1.cdna.all.fasta";

		CodePermutation P=new CodePermutation();
		P.loadDefaultcodeSet();
		new CodeEvaluation(P.calculateValues()).countBetterCodes(T);

		setWeightings(true, false, false, false);
		P=new CodePermutation();
		P.loadDefaultcodeSet();
		new CodeEvaluation(P.calculateValues()).countBetterCodes(T);

		setWeightings(false, true, false, false);
		P=new CodePermutation();
		P.loadDefaultcodeSet();
		new CodeEvaluation(P.calculateValues()).countBetterCodes(T);

		setWeightings(false, false, false, true);
		P=new CodePermutation();
		P.loadDefaultcodeSet();
		new CodeEvaluation(P.calculateValues()).countBetterCodes(T);

		setWeightings(true, true, false, false);
		P=new CodePermutation();
		P.loadDefaultcodeSet();
		new CodeEvaluation(P.calculateValues()).countBetterCodes(T);

		setWeightings(true, false, false, true);
		P=new CodePermutation();
		P.loadDefaultcodeSet();
		new CodeEvaluation(P.calculateValues()).countBetterCodes(T);

		setWeightings(false, true, false, true);
		P=new CodePermutation();
		P.loadDefaultcodeSet();
		new CodeEvaluation(P.calculateValues()).countBetterCodes(T);

		setWeightings(true, true, false, true);
		P=new CodePermutation();
		P.loadDefaultcodeSet();
		new CodeEvaluation(P.calculateValues()).countBetterCodes(T);


		DNASequence Seq1=conn.LoadFastaFile("Escherichia_coli.HUSEC2011CHR1.dna.chromosome.Chromosome");
		SequenceStats Stat2=new SequenceStats(Seq1.getSequenceAsString());
		baseAprioriWeights=Stat2.getBase_aPriori();
		tripletAprioriWeights=Stat2.getTriplet_aPriori();
		tripletTransitionWeights=Stat2.getTripletTransition();
		T="Tabelle 3.11 Reproduktion: Escherichia_coli.HUSEC2011CHR1.dna.chromosome.Chromosome";

        setWeightings(false, false, false, false);
		P=new CodePermutation();
		P.loadDefaultcodeSet();
		new CodeEvaluation(P.calculateValues()).countBetterCodes(T);

		setWeightings(true, false, false, false);
		P=new CodePermutation();
		P.loadDefaultcodeSet();
		new CodeEvaluation(P.calculateValues()).countBetterCodes(T);

		setWeightings(false, true, false, false);
		P=new CodePermutation();
		P.loadDefaultcodeSet();
		new CodeEvaluation(P.calculateValues()).countBetterCodes(T);

		setWeightings(false, false, false, true);
		P=new CodePermutation();
		P.loadDefaultcodeSet();
		new CodeEvaluation(P.calculateValues()).countBetterCodes(T);

		setWeightings(true, true, false, false);
		P=new CodePermutation();
		P.loadDefaultcodeSet();
		new CodeEvaluation(P.calculateValues()).countBetterCodes(T);

		setWeightings(true, false, false, true);
		P=new CodePermutation();
		P.loadDefaultcodeSet();
		new CodeEvaluation(P.calculateValues()).countBetterCodes(T);

		setWeightings(false, true, false, true);
		P=new CodePermutation();
		P.loadDefaultcodeSet();
		new CodeEvaluation(P.calculateValues()).countBetterCodes(T);

		setWeightings(true, true, false, true);
		P=new CodePermutation();
		P.loadDefaultcodeSet();
		new CodeEvaluation(P.calculateValues()).countBetterCodes(T);*/

		//ToDo: Gewichtung der Mutationen bei denen Stopcodons entstehen können in Codierenden Sequenzen und in nicht Codierenden Sequenzen
		//Entstehen bei Mutationen in der codierenden DNA mehr Stoppcodons als bei angenommener Gleichverteilung

		/*List<DNASequence> cDNA=conn.LoadMixedFileReadingframe("HomoSapiens_CCDS_Klaucke.fasta");
		System.out.println("Size: " +cDNA.size());
		SequenceStats_Coding Stat=new SequenceStats_Coding();
		for (DNASequence Seq : cDNA){
			Stat.ProcessSequence(Seq.getSequenceAsString());
		}
		Stat.FinalizeResults();
		baseAprioriWeights=Stat.getBase_aPriori();
		baseAprioriWeights=Stat.getBase_aPriori();
		tripletAprioriWeights=Stat.getTriplet_aPriori();
		baseTransitionWeights=Stat.getBaseTransition();
		tripletTransitionWeights=Stat.getTripletTransition();
		//None
		setWeightings(false,false,false,false);
		System.out.println("Unweighted Number of Possible Nonsense Mutations SNP: "+ Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_SNP()));
		System.out.println("Unweighted Number of Possible Nonsense Mutations Shift: "+ Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_Shift()));
		System.out.println("Sum: "+ToolMethods.getWeightedStopCodonFrequency_Overall());
		//Punktmutationen
		setWeightings(true,false,false,false);
		System.out.println("Weighted Number of Possible Nonsense Mutations SNP: "+ Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_SNP()));

		setWeightings(false,true,false,false);
		System.out.println("Weighted Number of Possible Nonsense Mutations SNP: "+ Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_SNP()));

		setWeightings(true,true,false,false);
		System.out.println("Weighted Number of Possible Nonsense Mutations SNP: "+ Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_SNP()));

		//Shift
		setWeightings(false,false,true,false);
		System.out.println("Weighted Number of Possible Nonsense Mutations Shift: "+ Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_Shift()));

		setWeightings(false,false,false,true);
		System.out.println("Weighted Number of Possible Nonsense Mutations Shift: "+ Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_Shift()));

		setWeightings(false,false,true,true);
		System.out.println("Weighted Number of Possible Nonsense Mutations Shift: "+ Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_Shift()));


		DNASequence Seq1=conn.LoadFastaFile("NC_000001.11");
		SequenceStats Stat2=new SequenceStats(Seq1.getSequenceAsString());
		baseAprioriWeights=Stat2.getBase_aPriori();
		tripletAprioriWeights=Stat2.getTriplet_aPriori();
		baseTransitionWeights=Stat2.getBaseTransition();
		tripletTransitionWeights=Stat2.getTripletTransition();

		setWeightings(true,false,false,false);
		System.out.println("Weighted Number of Possible Nonsense Mutations SNP: "+ Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_SNP()));

		setWeightings(false,true,false,false);
		System.out.println("Weighted Number of Possible Nonsense Mutations SNP: "+ Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_SNP()));

		setWeightings(true,true,false,false);
		System.out.println("Weighted Number of Possible Nonsense Mutations SNP: "+ Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_SNP()));

		//Shift
		setWeightings(false,false,true,false);
		System.out.println("Weighted Number of Possible Nonsense Mutations Shift: "+ Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_Shift()));

		setWeightings(false,false,false,true);
		System.out.println("Weighted Number of Possible Nonsense Mutations Shift: "+ Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_Shift()));

		setWeightings(false,false,true,true);
		System.out.println("Weighted Number of Possible Nonsense Mutations Shift: "+ Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_Shift()));*/

		//ToDo: WMS0 Vergleich abhängig von Transition/Transversion Bias, Tabelle 3.12

		//ToDo: Nukleotidverteilung in Codierenden Sequenzen vs ganzem Chromosom
		//ToDo: NT Gewichtung Differenz in Codierenden Sequenzen vs ganzem Chromosom
		/*List<DNASequence> cDNA=conn.LoadMixedFileReadingframe("HomoSapiens_CCDS_Klaucke.fasta");
		System.out.println("Size: " +cDNA.size());
		SequenceStats_Coding Stat=new SequenceStats_Coding();
		for (DNASequence Seq : cDNA){
			Stat.ProcessSequence(Seq.getSequenceAsString());
		}
		Stat.FinalizeResults();
		System.out.println("CCDS Matrizen");
		System.out.println("NA");
		ToolMethods.PrintMatrix(Stat.getBase_aPriori());
		System.out.println("TA");
		ToolMethods.PrintMatrix(Stat.getTriplet_aPriori());
		System.out.println("NT");
		ToolMethods.PrintMatrix(Stat.getBaseTransition());


		DNASequence Seq1=conn.LoadFastaFile("NC_000001.11");
		SequenceStats Stat2=new SequenceStats(Seq1.getSequenceAsString());
		System.out.println("Chromosom 1 Matrizen");
		System.out.println("NA");
		ToolMethods.PrintMatrix(Stat2.getBase_aPriori());
		System.out.println("TA");
		ToolMethods.PrintMatrix(Stat2.getTriplet_aPriori());
		System.out.println("NT");
		ToolMethods.PrintMatrix(Stat2.getBaseTransition());*/


		//ToDO: Histrgramm der Fehler mit Gewichtungen und ohne auf der Codierenden Sequenz
        //no Weightings
//		GeneCode code = new GeneCode();
//		setWeightings(false,false, false, false);
//		StabilityCalculator calc=new StabilityCalculator(code);
//		calc.setPrintHistogram(true);
//		calc.get_BaseDeviation(1);
//		calc.get_BaseDeviation(2);
//		calc.get_BaseDeviation(3);
//		calc.get_ShiftDeviation(1);
//		calc.get_ShiftDeviation(2);
//		//NA+TA+TT
//		List<DNASequence> cDNA=conn.LoadMixedFileReadingframe("HomoSapiens_CCDS_Klaucke.fasta");
//		System.out.println("Size: " +cDNA.size());
//		SequenceStats_Coding Stat=new SequenceStats_Coding();
//		for (DNASequence Seq : cDNA){
//			Stat.ProcessSequence(Seq.getSequenceAsString());
//		}
//		Stat.FinalizeResults();
//		baseAprioriWeights=Stat.getBase_aPriori();
//		baseAprioriWeights=Stat.getBase_aPriori();
//		tripletAprioriWeights=Stat.getTriplet_aPriori();
//		baseTransitionWeights=Stat.getBaseTransition();
//		tripletTransitionWeights=Stat.getTripletTransition();
//        setWeightings(true,true, false, true);
//        calc=new StabilityCalculator(code);
//        calc.setPrintHistogram(true);
//        calc.get_BaseDeviation(1);
//        calc.get_BaseDeviation(2);
//        calc.get_BaseDeviation(3);
//        calc.get_ShiftDeviation(1);
//        calc.get_ShiftDeviation(2);

		//ToDO: Messung der durchschnittlichen Länge der Sequenz nach einer Mutation bis zu einem Stoppcodon
		//Berechnung der Wahrscheinlichkeit nach jedem Triplett ein Stoppcodon zu erhalten
		//betrachtung der Länge bis maximal 20 Tripletts nach der Mutation
		//Basis: Triplettübergang zu Triplett (TT)
//		List<DNASequence> cDNA=conn.LoadMixedFile("HomoSapiens_CCDS_Klaucke.fasta");
//		System.out.println("Size: " +cDNA.size());
//		SequenceStats Stat=new SequenceStats();
//		for (DNASequence Seq : cDNA){
//			Stat.processSequence(Seq.getSequenceAsString());
//		}
//		Stat.finalizeResults();

		List<DNASequence> cDNA=conn.LoadMixedFileReadingframe("HomoSapiens_CCDS_Klaucke.fasta");
		System.out.println("Size: " +cDNA.size());
		SequenceStats_Coding Stat=new SequenceStats_Coding();
		for (DNASequence Seq : cDNA){
			Stat.processSequence(Seq.getSequenceAsString());
		}
		Stat.finalizeResults();

//		DNASequence Seq1=conn.LoadFastaFile("NC_000001.11");
//		SequenceStats Stat=new SequenceStats();
//		Stat.processSequence(Seq1.getSequenceAsString());
//		Stat.finalizeResults();

		StopCodonCalculator calc = new StopCodonCalculator(new GeneCode(),Stat.getBaseTransition(),Stat.getTriplet_aPriori());
		System.out.println(calc.getChainLengthAfterMuatation());




	/*	DNASequence Seq1=conn.LoadFastaFile("NC_000001.11");
		SequenceStats Stat2=new SequenceStats();
		Stat2.processSequence(Seq1.getSequenceAsString());
		Stat2.finalizeResults();
*/

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
		String Info="";
		if(ba)Info="[NA]";
		if(ta)Info=Info+"[TA]";
		if(bt)Info=Info+"[NT]";
		if(tt)Info=Info+"[TT]";
		if(Info.equals(""))Info="None";
		System.out.println("Weightings: "+Info);
	}

	private static void If(boolean ba) {
	}

}
