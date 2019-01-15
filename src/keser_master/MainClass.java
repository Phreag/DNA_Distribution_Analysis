package keser_master;

import java.util.Arrays;
import java.util.List;

import org.biojava.nbio.core.sequence.DNASequence;
import Objects.*;

public class MainClass {
    static GenBankConnection conn = new GenBankConnection();
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

    private static SequenceStatsCalculator loadSequenceStats(String GeneID, boolean readInReadingFrame) {
        DNASequence seq = conn.LoadFastaFile(GeneID);
        SequenceStatsCalculator stat = new SequenceStatsCalculator(readInReadingFrame, 0);
        stat.processSequence(seq.getSequenceAsString());
        return stat;
    }

    private static SequenceStatsCalculator loadSequenceStatsMixed(String sequenceFileIdentifier, boolean readInReadingFrame) {
        List<DNASequence> seqList;
        if (readInReadingFrame) {
            seqList = conn.LoadMixedFileReadingframe(sequenceFileIdentifier);
        } else {
            seqList = conn.LoadMixedFile(sequenceFileIdentifier);
        }
        System.out.println("Size: " + seqList.size());
        SequenceStatsCalculator stat = new SequenceStatsCalculator(readInReadingFrame, 0);
        for (DNASequence Seq : seqList) {
            stat.processSequence(Seq.getSequenceAsString());
        }
        return stat;
    }


    //This Program should be run with at least 7Gb of java heap space!
    public static void main(String[] args) {
        //compareNA_CCDS_CHR1();
        //stopCodonInspectionV1();
        //compareRandomCodesAcrossLifeforms();
        //stopCodonInspectionV2();
        stopCodonMarkovChainV1();



    }

    private static void compareNA_CCDS_CHR1() {
        /*
         * Nucleotidverteilungen auf Gesamter Sequenz vs Codierende Sequenz: Tabelle 3.2 (Claucke)
         */
        //CCDS
        SequenceStatsCalculator stat = loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true);
        System.out.println("CCDS NA Gewichte");
        ToolMethods.PrintMatrix(stat.getNucleotide_aPriori().getData());
        //Chromosom 1
        SequenceStatsCalculator stat2 = loadSequenceStats("NC_000001.11", false);
        System.out.println("Chromosom 1 NA Gewichte");
        ToolMethods.PrintMatrix(stat2.getNucleotide_aPriori().getData());
        /*Resultat
         * CCDS NA Gewichte T: 0.8668204368319248 C: 1.0654694621929026 A: 1.004639829034215 G: 1.0630702719409577
         * Chromosom 1 NA Gewichte T: 1.1670230561211625 C: 0.833995717634103 A: 1.164005262975987 G: 0.8349759632687476
         */
    }

    private static void stopCodonInspectionV1() {
        /*
         * Entstehen statistisch gesehen durch Gewichtungen mehr Stoppcondons?
         * Hypothese: Die DNA ist Optimiert dass es schnell zum Abbruch kommt
         * Vergleich: Komplette Sequenz vs Codierende Regionen
         */
        //Stoppcodons je Base bei Punktmutation Chromosom 1
        SequenceStatsCalculator stat = loadSequenceStats("NC_000001.11", false);
        WeightLoop loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getNucleotideTransition(), stat.gettripletTransition());
        System.out.println("Chromosom 1:");
        while (loop.moveNext()) {
            ToolMethods.getWeightedStopCodonFrequency_Overall();
        }
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
        SequenceStatsCalculator stat2 = loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true);
        WeightLoop loop2 = new WeightLoop(stat2.getNucleotide_aPriori(), stat2.getTriplet_aPriori(), stat2.getNucleotideTransition(), stat2.gettripletTransition());
        System.out.println("CCDS:");
        while (loop2.moveNext()) {
            ToolMethods.getWeightedStopCodonFrequency_Overall();
        }

    }

    private static void compareRandomCodesAcrossLifeforms() {
        //Tabelle 3.11 reproduzieren. Je 1 mal mit Homo Sapiens, E.Coli und Ciona intestinalis
        SequenceStatsCalculator stat = loadSequenceStatsMixed("Escherichia_coli.HUSEC2011CHR1.cdna.all.fasta", true);
        String T = "Tabelle 3.11 Reproduktion: Escherichia_coli.HUSEC2011CHR1.cdna.all.fasta";
        WeightLoop loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getNucleotideTransition(), stat.gettripletTransition());
        while (loop.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }
        SequenceStatsCalculator stat2 = loadSequenceStats("Escherichia_coli.HUSEC2011CHR1.dna.chromosome.Chromosome", false);
        T = "Tabelle 3.11 Reproduktion: Escherichia_coli.HUSEC2011CHR1.dna.chromosome.Chromosome";
        WeightLoop loop2 = new WeightLoop(stat2.getNucleotide_aPriori(), stat2.getTriplet_aPriori(), stat2.getNucleotideTransition(), stat2.gettripletTransition());
        while (loop2.moveNext()) {
            CodePermutation P = new CodePermutation();
            P.loadDefaultcodeSet();
            new CodeEvaluation(P.calculateValues()).countBetterCodes(T);
        }

    }

    private static void stopCodonInspectionV2() {
        //ToDo: Gewichtung der Mutationen bei denen Stopcodons entstehen können in Codierenden Sequenzen und in nicht Codierenden Sequenzen
        //Entstehen bei Mutationen in der codierenden DNA mehr Stoppcodons als bei angenommener Gleichverteilung

        SequenceStatsCalculator stat = loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true);
        WeightLoop loop = new WeightLoop(stat.getNucleotide_aPriori(), stat.getTriplet_aPriori(), stat.getNucleotideTransition(), stat.gettripletTransition());
        while (loop.moveNext()) {
            System.out.println("Number of Possible Nonsense Mutations SNP: " + Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_SNP()));
            System.out.println("Number of Possible Nonsense Mutations Shift: " + Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_Shift()));
            System.out.println("Sum: " + ToolMethods.getWeightedStopCodonFrequency_Overall());
        }

        SequenceStatsCalculator stat2 = loadSequenceStats("NC_000001.11", false);
        WeightLoop loop2 = new WeightLoop(stat2.getNucleotide_aPriori(), stat2.getTriplet_aPriori(), stat2.getNucleotideTransition(), stat2.gettripletTransition());
        while (loop2.moveNext()) {
            System.out.println("Number of Possible Nonsense Mutations SNP: " + Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_SNP()));
            System.out.println("Number of Possible Nonsense Mutations Shift: " + Arrays.toString(ToolMethods.getWeightedCountOfStopCodons_Shift()));
            System.out.println("Sum: " + ToolMethods.getWeightedStopCodonFrequency_Overall());
        }
    }

    private static void compareNT_CCDS_CHR1() {
        //ToDo: NT Gewichtung Differenz in Codierenden Sequenzen vs ganzem Chromosom
        SequenceStatsCalculator stat = loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true);
        System.out.println("CCDS Matrizen");
        System.out.println("NA");
        ToolMethods.PrintMatrix(stat.getNucleotide_aPriori().getData());
        System.out.println("TA");
        ToolMethods.PrintMatrix(stat.getTriplet_aPriori().getData());
        System.out.println("NT");
        ToolMethods.PrintMatrix(stat.getNucleotideTransition().getData());

        SequenceStatsCalculator stat2 = loadSequenceStats("NC_000001.11", false);
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
        SequenceStatsCalculator stat = loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true);
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
        SequenceStatsCalculator stat = loadSequenceStatsMixed("HomoSapiens_CCDS_Klaucke.fasta", true);
        MarkovChainForStopCodonsCalculator calc = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat.getNucleotideTransition(), stat.getTriplet_aPriori());
        System.out.println(calc.getChainLengthAfterMuatation());

        SequenceStatsCalculator stat2 = loadSequenceStats("NC_000001.11", false);
        MarkovChainForStopCodonsCalculator calc2 = new MarkovChainForStopCodonsCalculator(new GeneCode(), stat2.getNucleotideTransition(), stat2.getTriplet_aPriori());
        System.out.println(calc2.getChainLengthAfterMuatation());


    }
}
