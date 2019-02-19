package keser_master;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import Objects.*;
import org.biojava.nbio.core.sequence.DNASequence;

import javax.sound.midi.Sequence;

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


    //This Program should be run with at least 7Gb of java heap space!
    public static void main(String[] args) {
        //compareNA_CCDS_CHR1();
        //compareRandomCodesAcrossLifeforms();
        //nonsenseMutationCount();
        //stopCodonMarkovChainV2();
        //countStopCodonsInSequences();
        stopCodonMarkocChainCompareImpl();
        //getAverageCCDSSequenceLength();

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
        /*Resultat
         * CCDS NA Gewichte T: 0.8668204368319248 C: 1.0654694621929026 A: 1.004639829034215 G: 1.0630702719409577
         * Chromosom 1 NA Gewichte T: 1.1670230561211625 C: 0.833995717634103 A: 1.164005262975987 G: 0.8349759632687476
         */
    }

    private static void compareRandomCodesEColi() {
        //Tabelle 3.11 reproduzieren. Je 1 mal mit Homo Sapiens, E.Coli und Ciona intestinalis
        SequenceStatsCalculator stat = StatProvider.loadSequenceStatsMixed("Escherichia_coli.HUSEC2011CHR1.cdna.all.fasta", true);
        String T = "Tabelle 3.11 Reproduktion: Escherichia_coli.HUSEC2011CHR1.cdna.all.fasta";
        WeightLoop loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getNucleotideTransition(), stat.gettripletTransition());
        while (loop.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }
        SequenceStatsCalculator stat2 = StatProvider.loadSequenceStats("Escherichia_coli.HUSEC2011CHR1.dna.chromosome.Chromosome", false);
        T = "Tabelle 3.11 Reproduktion: Escherichia_coli.HUSEC2011CHR1.dna.chromosome.Chromosome";
        WeightLoop loop2 = new WeightLoop(stat2.getNucleotide_aPriori(), stat2.getTriplet_aPriori(), stat2.getNucleotideTransition(), stat2.gettripletTransition());
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
        WeightLoop loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getNucleotideTransition(), stat.gettripletTransition());
        while (loop.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }
        SequenceStatsCalculator stat2 = StatProvider.loadSequenceStats("NC_020166.2", false);
        T = "Tabelle 3.11 Reproduktion: Ciona_intestinalis Chromosom 1";
        WeightLoop loop2 = new WeightLoop(stat2.getNucleotide_aPriori(), stat2.getTriplet_aPriori(), stat2.getNucleotideTransition(), stat2.gettripletTransition());
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

        WeightLoop loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getNucleotideTransition(), stat.gettripletTransition());
        System.out.println("Nonsense Mutation Count on CCDS Homo Sapiens");
        while (loop.moveNext()) {
            System.out.println("#Nonsense Mutations SNP/Shift [Pos1, Pos2, Pos3][Left, Right]" + Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_SNP()) + Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_Shift()));
            System.out.println("Total sum of stopcodons: " + ToolMethods.getWeightedStopCodonFrequency_Overall());
        }
        WeightLoop loop2 = new WeightLoop(stat2.getNucleotide_aPriori(), stat2.getTriplet_aPriori(), stat2.getNucleotideTransition(), stat2.gettripletTransition());
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
        setWeightings(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), null, stat.gettripletTransition());
        calc = new StabilityCalculator(code);
        calc.setPrintHistogram(true);
        calc.get_BaseDeviation(1);
        calc.get_BaseDeviation(2);
        calc.get_BaseDeviation(3);
        calc.get_ShiftDeviation(1);
        calc.get_ShiftDeviation(2);
    }

    private static void stopCodonMarkovChainV1() {
        //ToDO: Messung der durchschnittlichen Länge der Sequenz nach einer Mutation bis zu einem Stoppcodon
        //Berechnung der Wahrscheinlichkeit nach jedem Triplett ein Stoppcodon zu erhalten
        //betrachtung der Länge bis maximal 20 Tripletts nach der Mutation
        //Basis: Triplettübergang zu Triplett (TT)
        SequenceStatsCalculator stat = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true);
        MarkovChainForStopCodonsCalculator calc = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat.getNucleotideTransition(), stat.getTriplet_aPriori());
        System.out.println(calc.getChainLengthAfterMuatation());

        SequenceStatsCalculator stat2 = StatProvider.loadSequenceStats("NC_000001.11", false);
        MarkovChainForStopCodonsCalculator calc2 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat2.getNucleotideTransition(), stat2.getTriplet_aPriori());
        System.out.println(calc2.getChainLengthAfterMuatation());


    }

    private static void stopCodonMarkovChainV2() {
        SequenceStatsCalculator stat0 = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true, 0);
        SequenceStatsCalculator stat1 = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true, 1);
        SequenceStatsCalculator stat2 = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true, 2);
        MarkovChainForStopCodonsCalculator calc0 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat0.getNucleotideTransition(), stat0.getTriplet_aPriori());
        System.out.println("CCDS Ofset 0: " + calc0.getChainLengthAfterMuatation());
        MarkovChainForStopCodonsCalculator calc1 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat1.getNucleotideTransition(), stat1.getTriplet_aPriori());
        System.out.println("CCDS Ofset 1: " + calc1.getChainLengthAfterMuatation());
        MarkovChainForStopCodonsCalculator calc2 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat2.getNucleotideTransition(), stat2.getTriplet_aPriori());
        System.out.println("CCDS Ofset 2: " + calc2.getChainLengthAfterMuatation());

        SequenceStatsCalculator statC0 = StatProvider.loadSequenceStats("NC_000001.11", true, 0);
        SequenceStatsCalculator statC1 = StatProvider.loadSequenceStats("NC_000001.11", true, 1);
        SequenceStatsCalculator statC2 = StatProvider.loadSequenceStats("NC_000001.11", true, 2);
        MarkovChainForStopCodonsCalculator calcC0 = new MarkovChainForStopCodonsCalculator(new GeneCode(), statC0.getNucleotideTransition(), statC0.getTriplet_aPriori());
        System.out.println("CCDS Ofset 0: " + calcC0.getChainLengthAfterMuatation());
        MarkovChainForStopCodonsCalculator calcC1 = new MarkovChainForStopCodonsCalculator(new GeneCode(), statC1.getNucleotideTransition(), statC1.getTriplet_aPriori());
        System.out.println("CCDS Ofset 1: " + calcC1.getChainLengthAfterMuatation());
        MarkovChainForStopCodonsCalculator calcC2 = new MarkovChainForStopCodonsCalculator(new GeneCode(), statC2.getNucleotideTransition(), statC2.getTriplet_aPriori());
        System.out.println("CCDS Ofset 2: " + calcC2.getChainLengthAfterMuatation());
    }

    private static void stopCodonMarkocChainCompareImpl() {
        //SequenceStatsCalculator stat0 = StatProvider.loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true, 0);
        SequenceStatsCalculator stat0 = StatProvider.loadSequenceStats("NC_000001.11", false, 0);
        TripletTransition TT = stat0.gettripletTransition();
        TripletApriori TA = stat0.getTriplet_aPriori();
        NucleotideTransition NT = stat0.getNucleotideTransition();
        GeneCode code = new GeneCode();

        List<Double> distancesRetDepth = new ArrayList<>();
        List<Double> distancesRetZero = new ArrayList<>();


        for (int depth = 0; depth < 200; depth++) {
            System.out.println("depth = "+depth);
            MarkovChainForStopCodonsCalculator calcDepth = new MarkovChainForStopCodonsCalculator(code, TT, TA, depth, false);
            MarkovChainForStopCodonsCalculator calcZero = new MarkovChainForStopCodonsCalculator(code, TT, TA, depth, true);
            distancesRetDepth.add(calcDepth.getChainLengthAfterMuatation());
            distancesRetZero.add(calcZero.getChainLengthAfterMuatation());
        }
        System.out.println("Distance return depth/zero");
        for (int depth = 0; depth < 200; depth++) {
            System.out.println(ToolMethods.df.format(distancesRetDepth.get(depth))+" / " +ToolMethods.df.format(distancesRetZero.get(depth)));
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

    private static void getAverageCCDSSequenceLength(){
        List<DNASequence> seqList = StatProvider.loadSequenceMixed("HomoSapiens_CCDS_Klaucke.fasta", true);
        int sum = 0;
        int count = seqList.size();
        for (DNASequence Seq : seqList) {
            String sequence = Seq.getSequenceAsString();
            int length=sequence.length();
            for (int i = 0 ; i < length-4; i = i + 3){
                if(Constants.isStopCodon(sequence.substring(i, i+3))){
                    System.out.println("Stop at Pos "+i+" / "+length);
                    count++;
                }
            }

            sum += length;

        }
        System.out.println("Average Sequence Length: "+(double)sum / (double)seqList.size());
        System.out.println("Average Sequence Length Adjusted: "+(double)sum / (double)count);
    }


}
