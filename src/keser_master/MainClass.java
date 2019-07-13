package keser_master;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import keser_master.MultiParam.CodeComparationMultiParam;
import keser_master.MultiParam.CodeEvaluationMultiParam;
import keser_master.MultiParam.CodeFinderMultiCharacteristics;
import keser_master.MultiParam.CodeFinderMultiScore;
import keser_master.Objects.*;
import org.biojava.nbio.core.sequence.DNASequence;

public class MainClass {
    public static NucleotideApriori nucleotideApriori;
    public static NucleotideTransition nucleotideTransition;
    public static TripletApriori tripletApriori;
    public static TripletTransition tripletTransition;
    public static int transitionTransversionBias = 1;

    public static void setWeightings(NucleotideApriori na, TripletApriori ta, NucleotideTransition nt, TripletTransition tt) {
        nucleotideApriori = na;
        tripletApriori = ta;
        nucleotideTransition = nt;
        tripletTransition = tt;
        System.out.println("Weightings: " + getConfigString());
    }

    public static String getConfigString() {
        String info = "";
        if (nucleotideApriori != null) info = "[NA]";
        if (tripletApriori != null) info = info + "[TA]";
        if (nucleotideTransition != null) info = info + "[NT]";
        if (tripletTransition != null) info = info + "[TT]";
        if (transitionTransversionBias != 1) {
            info = info + "[Bias = " + transitionTransversionBias + "]";
        }
        if (info.equals("")) info = "None";
        return info;
    }


    //This Program should be run with at least 8Gb of java heap space!
    //Stack trace needs to be increased for stop codon chain length calculation
    //used Arguments: -Xmx8G -Xss8m
    public static void main(String[] args) {
        //getTA_ZScores();
        //compareNA_CCDS_CHR1();
        //getTT2_ZScores();
        //compareRandomCodesAcrossLifeforms();
        //nonsenseMutationCount();
        //stopCodonMarkovChainV2();
        //countStopCodonsInSequences();
        //stopCodonMarkovChainCompareImpl();
        //getAverageCCDSSequenceLength();
        //generateRandomChromosome();
        //getAverageDistToEachCodon();
        //cleanTT2Weightings();
        //getAverageToEachCodonTA_Cleared();
        //millionHydropathyAndPolar();
        //millionHydropathyOnly();
        //millionCutOffHighDeltas();
        //millionOtherAminoAcidProperties();
        //millionMultiParam();
        //billionMultiParamTA_TT_Chr1();
        //runCodeFinderMultiCharacterstics();
        runCodeFinderMultiScore();

    }

    private static void getTA_ZScores() {
        SequenceStatsCalculator stat = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true, 0);
        System.out.println("Z-Scores");
        ToolMethods.PrintTripletTable(ToolMethods.calculateZScores(stat.getTriplet_aPriori().getData()), false);
        System.out.println("Weightings");
        ToolMethods.PrintTripletTable(stat.getTriplet_aPriori().getData(), false);


        stat = StatProvider.loadSequenceStats("NC_000001.11", false, 0);
        System.out.println("Z-Scores");
        ToolMethods.PrintTripletTable(ToolMethods.calculateZScores(stat.getTriplet_aPriori().getData()), false);
        System.out.println("Weightings");
        ToolMethods.PrintTripletTable(stat.getTriplet_aPriori().getData(), false);

    }

    private static void compareNA_CCDS_CHR1() {
        /*
         * Nucleotidverteilungen auf Gesamter Sequenz vs Codierende Sequenz: Tabelle 3.2 (Claucke)
         */
        //CCDS
        SequenceStatsCalculator stat = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true);
        System.out.println("CCDS NA Gewichte");
        ToolMethods.PrintMatrix(stat.getNucleotide_aPriori().getData());
        //Chromosom 1
        SequenceStatsCalculator stat2 = StatProvider.loadSequenceStats("NC_000001.11", false);
        System.out.println("Chromosom 1 NA Gewichte");
        ToolMethods.PrintMatrix(stat2.getNucleotide_aPriori().getData());
    }

    private static void getTT2_ZScores(){
        SequenceStatsCalculator stat = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true);
        System.out.println("CCDS NA Gewichte");
        ToolMethods.printTT2Table(ToolMethods.calculateZScores(stat.getTripletTransition().getDataTriplet()));
        //Chromosom 1
        SequenceStatsCalculator stat2 = StatProvider.loadSequenceStats("NC_000001.11", false);
        System.out.println("Chromosom 1 NA Gewichte");
        ToolMethods.printTT2Table(ToolMethods.calculateZScores(stat2.getTripletTransition().getDataTriplet()));
    }

    private static void compareRandomCodesEColi() {
        //Tabelle 3.11 reproduzieren. Je 1 mal mit Homo Sapiens, E.Coli und Ciona intestinalis
        SequenceStatsCalculator stat = StatProvider.loadSequenceStatsMixed("Escherichia_coli.HUSEC2011CHR1.cdna.all.fasta", true);
        String T = "Tabelle 3.11 Reproduktion: Escherichia_coli.HUSEC2011CHR1.cdna.all.fasta";
        WeightLoop loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getNucleotideTransition(), stat.getTripletTransition());
        while (loop.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }
        SequenceStatsCalculator stat2 = StatProvider.loadSequenceStats("Escherichia_coli.HUSEC2011CHR1.dna.chromosome.Chromosome", false);
        T = "Tabelle 3.11 Reproduktion: Escherichia_coli.HUSEC2011CHR1.dna.chromosome.Chromosome";
        WeightLoop loop2 = new WeightLoop(stat2.getNucleotide_aPriori(), stat2.getTriplet_aPriori(), stat2.getNucleotideTransition(), stat2.getTripletTransition());
        while (loop2.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }
    }

    private static void compareRandomCodesCionaI() {
        //Tabelle 3.11 reproduzieren. Je 1 mal mit Homo Sapiens, E.Coli und Ciona intestinalis
        SequenceStatsCalculator stat = StatProvider.loadSequenceStatsMixed("Ciona_intestinalis_CCDS.fasta", true);
        String T = "Tabelle 3.11 Reproduktion: Ciona_intestinalis_CCDS.fasta";
        WeightLoop loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getNucleotideTransition(), stat.getTripletTransition());
        while (loop.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }
        SequenceStatsCalculator stat2 = StatProvider.loadSequenceStats("NC_020166.2", false);
        T = "Tabelle 3.11 Reproduktion: Ciona_intestinalis Chromosom 1";
        WeightLoop loop2 = new WeightLoop(stat2.getNucleotide_aPriori(), stat2.getTriplet_aPriori(), stat2.getNucleotideTransition(), stat2.getTripletTransition());
        while (loop2.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }
    }

    private static void nonsenseMutationCount() {
        //ToDo: Gewichtung der Mutationen bei denen Stopcodons entstehen können in Codierenden Sequenzen und in nicht Codierenden Sequenzen
        //Entstehen bei Mutationen in der codierenden DNA mehr Stoppcodons als bei angenommener Gleichverteilung

        SequenceStatsCalculator stat = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true);
        SequenceStatsCalculator stat2 = StatProvider.loadSequenceStats("NC_000001.11", false);

        WeightLoop loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getNucleotideTransition(), stat.getTripletTransition());
        System.out.println("Nonsense Mutation Count on CCDS Homo Sapiens");
        while (loop.moveNext()) {
            System.out.println("#Nonsense Mutations SNP/Shift [Pos1, Pos2, Pos3][Left, Right]" + Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_SNP()) + Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_Shift()));
            System.out.println("Total sum of stopcodons: " + ToolMethods.getWeightedStopCodonFrequency_Overall());
        }
        WeightLoop loop2 = new WeightLoop(stat2.getNucleotide_aPriori(), stat2.getTriplet_aPriori(), stat2.getNucleotideTransition(), stat2.getTripletTransition());
        System.out.println("Nonsense Mutation Count on Chromosome 1 Homo Sapiens");
        while (loop2.moveNext()) {
            System.out.println("#Nonsense Mutations SNP/Shift [Pos1, Pos2, Pos3][Left, Right]" + Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_SNP()) + Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_Shift()));
            System.out.println("Total sum of stopcodons: " + ToolMethods.getWeightedStopCodonFrequency_Overall());
        }
    }

    private static void compareNT_CCDS_CHR1() {
        //ToDo: NT Gewichtung Differenz in Codierenden Sequenzen vs ganzem Chromosom
        SequenceStatsCalculator stat = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true);
        System.out.println("CCDS Matrizen");
        System.out.println("NA");
        ToolMethods.PrintMatrix(stat.getNucleotide_aPriori().getData());
        System.out.println("TA");
        ToolMethods.PrintMatrix(stat.getTriplet_aPriori().getData());
        System.out.println("NT");
        ToolMethods.PrintMatrix(stat.getNucleotideTransition().getData());

        SequenceStatsCalculator stat2 = StatProvider.loadSequenceStats("NC_000001.11", false);
        System.out.println("Chromosom 1 Matrizen");
        System.out.println("NA");
        ToolMethods.PrintMatrix(stat2.getNucleotide_aPriori().getData());
        System.out.println("TA");
        ToolMethods.PrintMatrix(stat2.getTriplet_aPriori().getData());
        System.out.println("NT");
        ToolMethods.PrintMatrix(stat2.getNucleotideTransition().getData());
    }

    private static void calculateErrorHistogram() {
        //ToDO: Histrgramm der Fehler mit Gewichtungen und ohne auf der Codierenden Sequenz
        //no Weightings
        GeneCode code = new GeneCode();
        StabilityCalculator calc = new StabilityCalculator(code);
        calc.setPrintHistogram(true);
        calc.get_BaseDeviation(1);
        calc.get_BaseDeviation(2);
        calc.get_BaseDeviation(3);
        calc.get_ShiftDeviation(1);
        calc.get_ShiftDeviation(2);
        //NA+TA+TT
        SequenceStatsCalculator stat = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true);
        setWeightings(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), null, stat.getTripletTransition());
        calc = new StabilityCalculator(code);
        calc.setPrintHistogram(true);
        calc.get_BaseDeviation(1);
        calc.get_BaseDeviation(2);
        calc.get_BaseDeviation(3);
        calc.get_ShiftDeviation(1);
        calc.get_ShiftDeviation(2);
    }

    private static void stopCodonMarkovChain() {
        SequenceStatsCalculator stat0 = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true, 0);
        SequenceStatsCalculator stat1 = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true, 1);
        SequenceStatsCalculator stat2 = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true, 2);
        MarkovChainForStopCodonsCalculator calcTT0 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat0.getTripletTransition(), stat0.getTriplet_aPriori(), 20000, false);
        MarkovChainForStopCodonsCalculator calcTT1 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat1.getTripletTransition(), stat1.getTriplet_aPriori(), 20000, false);
        MarkovChainForStopCodonsCalculator calcTT2 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat2.getTripletTransition(), stat2.getTriplet_aPriori(), 20000, false);
        MarkovChainForStopCodonsCalculator calcNT0 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat0.getNucleotideTransition(), stat0.getTriplet_aPriori(), 20000, false);
        MarkovChainForStopCodonsCalculator calcNT1 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat1.getNucleotideTransition(), stat1.getTriplet_aPriori(), 20000, false);
        MarkovChainForStopCodonsCalculator calcNT2 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat2.getNucleotideTransition(), stat2.getTriplet_aPriori(), 20000, false);
        System.out.println("CCDS TT2 Offset 0: " + calcTT0.getChainLengthAfterMuatation());
        System.out.println("CCDS TT2 Offset 1: " + calcTT1.getChainLengthAfterMuatation());
        System.out.println("CCDS TT2 Offset 2: " + calcTT2.getChainLengthAfterMuatation());
        System.out.println("CCDS NT Offset 0: " + calcNT0.getChainLengthAfterMuatation());
        System.out.println("CCDS NT Offset 1: " + calcNT1.getChainLengthAfterMuatation());
        System.out.println("CCDS NT Offset 2: " + calcNT2.getChainLengthAfterMuatation());

        stat0 = StatProvider.loadSequenceStats("RandomSequence", true, 0);
        stat1 = StatProvider.loadSequenceStats("RandomSequence", true, 1);
        stat2 = StatProvider.loadSequenceStats("RandomSequence", true, 2);
        calcTT0 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat0.getTripletTransition(), stat0.getTriplet_aPriori(), 20000, false);
        calcTT1 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat1.getTripletTransition(), stat1.getTriplet_aPriori(), 20000, false);
        calcTT2 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat2.getTripletTransition(), stat2.getTriplet_aPriori(), 20000, false);
        calcNT0 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat0.getNucleotideTransition(), stat0.getTriplet_aPriori(), 20000, false);
        calcNT1 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat1.getNucleotideTransition(), stat1.getTriplet_aPriori(), 20000, false);
        calcNT2 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat2.getNucleotideTransition(), stat2.getTriplet_aPriori(), 20000, false);
        System.out.println("Random TT2 Offset 0: " + calcTT0.getChainLengthAfterMuatation());
        System.out.println("Random TT2 Offset 1: " + calcTT1.getChainLengthAfterMuatation());
        System.out.println("Random TT2 Offset 2: " + calcTT2.getChainLengthAfterMuatation());
        System.out.println("Random NT Offset 0: " + calcNT0.getChainLengthAfterMuatation());
        System.out.println("Random NT Offset 1: " + calcNT1.getChainLengthAfterMuatation());
        System.out.println("Random NT Offset 2: " + calcNT2.getChainLengthAfterMuatation());

        stat0 = StatProvider.loadSequenceStats("NC_000001.11", true, 0);
        stat1 = StatProvider.loadSequenceStats("NC_000001.11", true, 1);
        stat2 = StatProvider.loadSequenceStats("NC_000001.11", true, 2);
        calcTT0 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat0.getTripletTransition(), stat0.getTriplet_aPriori(), 20000, false);
        calcTT1 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat1.getTripletTransition(), stat1.getTriplet_aPriori(), 20000, false);
        calcTT2 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat2.getTripletTransition(), stat2.getTriplet_aPriori(), 20000, false);
        calcNT0 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat0.getNucleotideTransition(), stat0.getTriplet_aPriori(), 20000, false);
        calcNT1 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat1.getNucleotideTransition(), stat1.getTriplet_aPriori(), 20000, false);
        calcNT2 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat2.getNucleotideTransition(), stat2.getTriplet_aPriori(), 20000, false);
        System.out.println("CHR1 TT2 Offset 0: " + calcTT0.getChainLengthAfterMuatation());
        System.out.println("CHR1 TT2 Offset 1: " + calcTT1.getChainLengthAfterMuatation());
        System.out.println("CHR1 TT2 Offset 2: " + calcTT2.getChainLengthAfterMuatation());
        System.out.println("CHR1 NT Offset 0: " + calcNT0.getChainLengthAfterMuatation());
        System.out.println("CHR1 NT Offset 1: " + calcNT1.getChainLengthAfterMuatation());
        System.out.println("CHR1 NT Offset 2: " + calcNT2.getChainLengthAfterMuatation());

    }

    private static void stopCodonMarkovChainCompareImpl() {
        try {
            SequenceStatsCalculator stat0 = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true, 0);
            TripletTransition TT = stat0.getTripletTransition();
            TripletApriori TA = stat0.getTriplet_aPriori();
            FileWriter fw = new FileWriter(new File("data/CompareTreeImpl.csv"), false);
            GeneCode code = new GeneCode();


            fw.write("Search Distance; return Depth; return Zero\n");
            for (int depth = 0; depth < 5000; depth = depth + 50) {
                MarkovChainForStopCodonsCalculator calcDepth = new MarkovChainForStopCodonsCalculator(code, TT, TA, depth, false);
                MarkovChainForStopCodonsCalculator calcZero = new MarkovChainForStopCodonsCalculator(code, TT, TA, depth, true);
                String line = depth + "; " + ToolMethods.df.format(calcDepth.getChainLengthAfterMuatation()) + "; " + ToolMethods.df.format(calcZero.getChainLengthAfterMuatation());
                System.out.println(line);
                fw.write(line + "\n");
            }
            fw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void countStopCodonsInSequences() {
        DNASequence chr1 = StatProvider.loadSequence("NC_000001.11");
        List<DNASequence> ccds = StatProvider.loadSequenceMixed("HomoSapiens_CCDS_Klaucke.fasta", true);
        System.out.println("Chromosome1 Stopcodonfrequency Offset 0: " + ToolMethods.df_long.format(ToolMethods.getStopCodonFrequency(Collections.singletonList(chr1), true, 0)));
        System.out.println("Chromosome1 Stopcodonfrequency Offset 1: " + ToolMethods.df_long.format(ToolMethods.getStopCodonFrequency(Collections.singletonList(chr1), true, 1)));
        System.out.println("Chromosome1 Stopcodonfrequency Offset 2: " + ToolMethods.df_long.format(ToolMethods.getStopCodonFrequency(Collections.singletonList(chr1), true, 2)));
        System.out.println("Chromosome1 Stopcodonfrequency Non Readingframe: " + ToolMethods.df_long.format(ToolMethods.getStopCodonFrequency(Collections.singletonList(chr1), false, 0)));

        System.out.println("CCDS Stopcodonfrequency Offset 0: " + ToolMethods.df_long.format(ToolMethods.getStopCodonFrequency(ccds, true, 0)));
        System.out.println("CCDS Stopcodonfrequency Offset 1: " + ToolMethods.df_long.format(ToolMethods.getStopCodonFrequency(ccds, true, 1)));
        System.out.println("CCDS Stopcodonfrequency Offset 2: " + ToolMethods.df_long.format(ToolMethods.getStopCodonFrequency(ccds, true, 2)));
        System.out.println("CCDS Stopcodonfrequency Non Readingframe: " + ToolMethods.df_long.format(ToolMethods.getStopCodonFrequency(ccds, false, 0)));
    }

    private static void getAverageCCDSSequenceLength() {
        List<DNASequence> seqList = StatProvider.loadSequenceMixed("HomoSapiens_CCDS_Klaucke.fasta", true);
        int sum = 0;
        int count = 0;
        for (DNASequence Seq : seqList) {
            String sequence = Seq.getSequenceAsString();
            int length = sequence.length();
            for (int i = 0; i < length - 4; i = i + 3) {
                if (Constants.isStopCodon(sequence.substring(i, i + 3))) {
                    System.out.println("Stop at Pos " + i + " / " + length);
                    count++;
                }
            }

            sum += length;

        }
        System.out.println("Stopcodons in Readingframe: " + count);
        System.out.println("Average Sequence Length: " + (double) sum / (double) seqList.size());
        System.out.println("Average Sequence Length Adjusted: " + (double) sum / (double) (count + seqList.size()));
    }

    private static void generateRandomChromosome() {
        try {
            FileWriter fw = new FileWriter(new File("data/RandomSequence.fasta"), false);
            BufferedWriter bw = new BufferedWriter(fw);
            XORShift_Random rnd = new XORShift_Random();
            bw.write(">gi|0000000|ref|RandomSequence| random sequence generated with Xorshift\n");
            int[] nucCount = new int[4];
            for (int i = 0; i < 250000000; i++) {
                int nucleotide = rnd.getInt(4);
                nucCount[nucleotide]++;
                bw.write(Constants.Bases[nucleotide]);
                if ((i + 1) % 70 == 0) {
                    bw.write("\n");
                }
            }
            System.out.println(nucCount[0] + " / " + nucCount[1] + " / " + nucCount[2] + " / " + nucCount[3]);
            bw.close();
        } catch (IOException e) {
            System.out.println("Filewriter Error");
            e.printStackTrace();
        }
    }

    private static void getAverageDistToEachCodon() {
        SequenceStatsCalculator stat0 = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true, 0);
        TripletTransition TT = stat0.getTripletTransition();
        TripletApriori TA = stat0.getTriplet_aPriori();
        GeneCode code = new GeneCode();

        double[][][] Values = new double[4][4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    String codon = Constants.Bases[i] + Constants.Bases[j] + Constants.Bases[k];
                    MarkovChainForStopCodonsCalculator calc = new MarkovChainForStopCodonsCalculator(code, TT, TA, 20000, codon);
                    Values[i][j][k] = calc.getChainLengthAfterMuatation();
                }
            }
        }
        System.out.println("Average Chain length to specific codon (CCDS):");
        ToolMethods.PrintTripletTable(Values, true);
        ToolMethods.PrintTripletTable(ToolMethods.calculateZScores(Values),true);

        SequenceStatsCalculator statC0 = StatProvider.loadSequenceStats("NC_000001.11", false, 0);
        TT = statC0.getTripletTransition();
        TA = statC0.getTriplet_aPriori();

        Values = new double[4][4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    String codon = Constants.Bases[i] + Constants.Bases[j] + Constants.Bases[k];
                    MarkovChainForStopCodonsCalculator calc = new MarkovChainForStopCodonsCalculator(code, TT, TA, 20000, codon);
                    Values[i][j][k] = calc.getChainLengthAfterMuatation();
                }
            }
        }
        System.out.println("Average Chain length to specific codon (Chr1):");
        ToolMethods.PrintTripletTable(Values, true);
        ToolMethods.PrintTripletTable(ToolMethods.calculateZScores(Values),true);

        SequenceStatsCalculator statR0 = StatProvider.loadSequenceStats("RandomSequence", false, 0);
        TT = statR0.getTripletTransition();
        TA = statR0.getTriplet_aPriori();
        Values = new double[4][4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    String codon = Constants.Bases[i] + Constants.Bases[j] + Constants.Bases[k];
                    MarkovChainForStopCodonsCalculator calc = new MarkovChainForStopCodonsCalculator(code, TT, TA, 20000, codon);
                    Values[i][j][k] = calc.getChainLengthAfterMuatation();
                }
            }
        }
        System.out.println("Average Chain length to specific codon (Random):");
        ToolMethods.PrintTripletTable(Values, true);
        ToolMethods.PrintTripletTable(ToolMethods.calculateZScores(Values),true);
    }

    private static void getAverageToEachCodonTA_Cleared() {
        SequenceStatsCalculator stat0 = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true, 0);
        TripletTransition TT = stat0.getTripletTransition();
        TripletApriori TA = stat0.getTriplet_aPriori();
        TT.cleanTriplettApriori(TA);
        GeneCode code = new GeneCode();

        double[][][] Values = new double[4][4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    String codon = Constants.Bases[i] + Constants.Bases[j] + Constants.Bases[k];
                    MarkovChainForStopCodonsCalculator calc = new MarkovChainForStopCodonsCalculator(code, TT, TA, 20000, codon);
                    Values[i][j][k] = calc.getChainLengthAfterMuatation();
                }
            }
        }
        System.out.println("Average Chain length to specific codon (CCDS & Cleared):");
        ToolMethods.PrintTripletTable(Values, true);
        System.out.println("Z-Scores");
        ToolMethods.PrintTripletTable(ToolMethods.calculateZScores(Values), true);

        SequenceStatsCalculator statC0 = StatProvider.loadSequenceStats("NC_000001.11", false, 0);
        TT = statC0.getTripletTransition();
        TA = statC0.getTriplet_aPriori();
        TT.cleanTriplettApriori(TA);

        Values = new double[4][4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    String codon = Constants.Bases[i] + Constants.Bases[j] + Constants.Bases[k];
                    MarkovChainForStopCodonsCalculator calc = new MarkovChainForStopCodonsCalculator(code, TT, TA, 20000, codon);
                    Values[i][j][k] = calc.getChainLengthAfterMuatation();
                }
            }
        }
        System.out.println("Average Chain length to specific codon (Chr1 & Cleared):");
        ToolMethods.PrintTripletTable(Values, true);
        System.out.println("Z-Scores");
        ToolMethods.PrintTripletTable(ToolMethods.calculateZScores(Values), true);

        SequenceStatsCalculator statR0 = StatProvider.loadSequenceStats("RandomSequence", false, 0);
        TT = statR0.getTripletTransition();
        TA = statR0.getTriplet_aPriori();
        TT.cleanTriplettApriori(TA);
        Values = new double[4][4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    String codon = Constants.Bases[i] + Constants.Bases[j] + Constants.Bases[k];
                    MarkovChainForStopCodonsCalculator calc = new MarkovChainForStopCodonsCalculator(code, TT, TA, 20000, codon);
                    Values[i][j][k] = calc.getChainLengthAfterMuatation();
                }
            }
        }
        System.out.println("Average Chain length to specific codon (Random & Cleared):");
        ToolMethods.PrintTripletTable(Values, true);
        System.out.println("Z-Scores");
        ToolMethods.PrintTripletTable(ToolMethods.calculateZScores(Values), true);
    }

    private static void cleanTT2Weightings() {
        SequenceStatsCalculator stat0 = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true, 0);
        TripletTransition TT = stat0.getTripletTransition();
        TripletApriori TA = stat0.getTriplet_aPriori();
        TT.cleanTriplettApriori(TA);
        GeneCode code = new GeneCode();
        double[][][] Values = new double[4][4][4];
        double[][][] Values2 = new double[4][4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    String codon = Constants.Bases[i] + Constants.Bases[j] + Constants.Bases[k];
                    MarkovChainForStopCodonsCalculator calc = new MarkovChainForStopCodonsCalculator(code, TT, TA, 20000, codon);
                    Values[i][j][k] = calc.getChainLengthAfterMuatation();
                    Values2[i][j][k] = calc.getChainLengthAfterMuatation();
                }
            }
        }
        System.out.println("Average Chain length to specific codon (CCDS):");
        ToolMethods.PrintTripletTable(Values, true);
        System.out.println("Z-Scores (CCDS):");
        ToolMethods.PrintTripletTable(ToolMethods.calculateZScores(Values2), true);

    }

    private static void millionHydropathyAndPolar() {
        //ref
        CodePermutation Pe = new CodePermutation();
        Pe.loadDefaultcodeSet();
        new CodeEvaluation(Pe.calculateValues()).countBetterCodes("No Weightings, Default settings");

        //enable hydropathy
        Constants.hydropathyEnabled = true;
        Constants.normalizeBySigma = true;
        SequenceStatsCalculator stat = StatProvider.loadSequenceStats("NC_000001.11", false, 0);
        String T = "Chromosome 1 hydropathy + polar (normalized by Sigma)";
        WeightLoop loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getNucleotideTransition(), stat.getTripletTransition());
        while (loop.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }

        stat = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true, 0);
        T = "CCDS hydropathy + polar (normalized by Sigma)";
        loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getNucleotideTransition(), stat.getTripletTransition());
        while (loop.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }
    }

    private static void millionHydropathyOnly() {
        //enable hydropathy
        Constants.hydropathyEnabled = true;
        Constants.polarReqEnabled = false;
        SequenceStatsCalculator stat = StatProvider.loadSequenceStats("NC_000001.11", false, 0);
        String T = "Chromosome 1 polar";
        WeightLoop loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getTripletTransition());
        while (loop.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }

        stat = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true, 0);
        T = "CCDS polar";
        loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getTripletTransition());
        while (loop.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }
    }

    private static void millionCutOffHighDeltas(){
        Constants.cutOff = 66;
        SequenceStatsCalculator stat = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true, 0);
        String T = "CCDS polar, cutOff: 66";
        WeightLoop loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getTripletTransition());
        while (loop.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }
    }

    private static void millionOtherAminoAcidProperties(){
        Constants.polarReqEnabled = false;
        //molecular Weight
        Constants.molVolEnabled = true;
        SequenceStatsCalculator stat = StatProvider.loadSequenceStats("NC_000001.11", false, 0);
        String T = "Chromosome 1 movlolume";
        WeightLoop loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getTripletTransition());
        while (loop.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }
        stat = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true, 0);
        T = "CCDS molVolume";
        loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getTripletTransition());
        while (loop.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }

        //molecular Weight
        Constants.molVolEnabled = false;
        Constants.molWeightEnabled = true;
        stat = StatProvider.loadSequenceStats("NC_000001.11", false, 0);
        T = "Chromosome 1 molWeight";
        loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getTripletTransition());
        while (loop.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }
        stat = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true, 0);
        T = "CCDS molWeight";
        loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getTripletTransition());
        while (loop.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }

        //pKa
        Constants.molWeightEnabled=false;
        Constants.pKaEnabled=true;
        stat = StatProvider.loadSequenceStats("NC_000001.11", false, 0);
        T = "Chromosome 1 pKa";
        loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getTripletTransition());
        while (loop.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }
        stat = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true, 0);
        T = "CCDS pKa";
        loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getTripletTransition());
        while (loop.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }

        //pKb
        Constants.pKaEnabled=false;
        Constants.pKbEnabled=true;
        stat = StatProvider.loadSequenceStats("NC_000001.11", false, 0);
        T = "Chromosome 1 pKb";
        loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getTripletTransition());
        while (loop.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }
        stat = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true, 0);
        T = "CCDS pKb";
        loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getTripletTransition());
        while (loop.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }

        //pI
        Constants.pKbEnabled=false;
        Constants.pIEnabled=true;
        stat = StatProvider.loadSequenceStats("NC_000001.11", false, 0);
        T = "Chromosome 1 pI";
        loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getTripletTransition());
        while (loop.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }
        stat = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true, 0);
        T = "CCDS pI";
        loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getTripletTransition());
        while (loop.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }
    }

    private static void millionMultiParam(){
        //generate Codes
        CodeComparationMultiParam comp = new CodeComparationMultiParam();
        for (int i = 0 ; i < 1000 ; i++){
            comp.generateCodeSet("RandomCodes_"+i);
        }

        SequenceStatsCalculator stat = StatProvider.loadSequenceStats("NC_000001.11", false, 0);
        System.gc();
        WeightLoop loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getTripletTransition());
        while (loop.moveNext()) {
            calculateMillionMultiParam("Chromosome 1");
        }

        stat = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true, 0);
        System.gc();
        loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getTripletTransition());
        while (loop.moveNext()) {
            calculateMillionMultiParam("CCDS");
        }
    }

    private static void calculateMillionMultiParam (String title){
        CodeComparationMultiParam comp = new CodeComparationMultiParam();
        CodeEvaluationMultiParam evaluation = new CodeEvaluationMultiParam(null);
//        comp.loadDefaultcodeSet();
//        evaluation.setValues(comp.calculateValues());
//        evaluation.countBetterCodes(title);
        for (int i = 0 ; i < 1000 ; i++){
            comp.loadCodeSet("RandomCodes_"+i);
            evaluation.setValues(comp.calculateValues());
            evaluation.countBetterCodes(title);
        }
        evaluation.printSummary();
    }

    private static void billionMultiParamTA_TT_Chr1(){
        CodeComparationMultiParam comp = new CodeComparationMultiParam();
        for (int i = 0 ; i < 10000 ; i++){
            comp.generateCodeSet("RandomCodes_"+i);
        }

        SequenceStatsCalculator stat = StatProvider.loadSequenceStats("NC_000001.11", false, 0);
        System.gc();
        setWeightings(null,stat.getTriplet_aPriori(),null,stat.getTripletTransition());
        comp = new CodeComparationMultiParam();
        CodeEvaluationMultiParam evaluation = new CodeEvaluationMultiParam(null);
        for (int i = 0 ; i < 10000 ; i++){
            comp.loadCodeSet("RandomCodes_"+i);
            evaluation.setValues(comp.calculateValues());
            evaluation.countBetterCodes("Chromosome 1");
        }
        evaluation.printSummary();
    }

    private static void runCodeFinderMultiCharacterstics(){
        SequenceStatsCalculator stat = StatProvider.loadSequenceStats("NC_000001.11", false, 0);
        System.gc();
        setWeightings(null,stat.getTriplet_aPriori(),null,stat.getTripletTransition());
        CodeComparationMultiParam.loadDefaultcodeSet();
        CodeFinderMultiCharacteristics greedy = new CodeFinderMultiCharacteristics();
        greedy.RunCodeFinder(50);
    }

    private static void runCodeFinderMultiScore(){
        SequenceStatsCalculator stat = StatProvider.loadSequenceStats("NC_000001.11", false, 0);
        System.gc();
        setWeightings(null,stat.getTriplet_aPriori(),null,stat.getTripletTransition());
        CodeComparationMultiParam.loadDefaultcodeSet();
        CodeFinderMultiScore greedy = new CodeFinderMultiScore();
        greedy.RunCodeFinder(50);
    }
}
