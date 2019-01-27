package keser_master;

import Objects.SequenceStatsCalculator;
import org.biojava.nbio.core.sequence.DNASequence;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class StatProvider {
    static GenBankConnection conn = new GenBankConnection();
    private static Map<String, DNASequence> singleSequenceCache = new HashMap<>();
    private static Map<String, List<DNASequence>> mixedSequenceCache = new HashMap<>();

    public static SequenceStatsCalculator loadSequenceStats(String GeneID, boolean readInReadingFrame) {
        return loadSequenceStats(GeneID, readInReadingFrame, 0);
    }

    public static SequenceStatsCalculator loadSequenceStats(String GeneID, boolean readInReadingFrame, int offset) {
        DNASequence seq = loadSequence(GeneID);
        SequenceStatsCalculator stat = new SequenceStatsCalculator(readInReadingFrame, offset);
        stat.processSequence(seq.getSequenceAsString());
        return stat;
    }

    public static SequenceStatsCalculator loadSequenceStatsMixed(String sequenceFileIdentifier, boolean readInReadingFrame) {
        return loadSequenceStatsMixed(sequenceFileIdentifier, readInReadingFrame, 0);
    }

    public static SequenceStatsCalculator loadSequenceStatsMixed(String sequenceFileIdentifier, boolean readInReadingFrame, int offset) {
        List<DNASequence> seqList = loadSequenceMixed(sequenceFileIdentifier, readInReadingFrame);
        System.out.println("Size: " + seqList.size());
        SequenceStatsCalculator stat = new SequenceStatsCalculator(readInReadingFrame, offset);
        for (DNASequence Seq : seqList) {
            stat.processSequence(Seq.getSequenceAsString());
        }
        return stat;
    }

    public static List<DNASequence> loadSequenceMixed (String sequenceFileIdentifier, boolean readInReadingFrame){
        if (readInReadingFrame) {
            if (!mixedSequenceCache.containsKey(sequenceFileIdentifier + "_ReadingFrame")) {
                mixedSequenceCache.put(sequenceFileIdentifier + "_ReadingFrame", conn.LoadMixedFileReadingframe(sequenceFileIdentifier));
            }
            return mixedSequenceCache.get(sequenceFileIdentifier + "_ReadingFrame");
        } else {
            if (!mixedSequenceCache.containsKey(sequenceFileIdentifier)) {
                mixedSequenceCache.put(sequenceFileIdentifier, conn.LoadMixedFile(sequenceFileIdentifier));
            }
            return mixedSequenceCache.get(sequenceFileIdentifier);
        }
    }

    public static DNASequence loadSequence (String GeneID){
        if (!singleSequenceCache.containsKey(GeneID)) {
            singleSequenceCache.put(GeneID, conn.LoadFastaFile(GeneID));
        }
        return singleSequenceCache.get(GeneID);
    }
}
