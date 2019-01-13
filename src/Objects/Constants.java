package Objects;

import java.security.InvalidParameterException;

public class Constants {

    public static String[] Bases = {"T", "C", "A", "G"};

    //Holds data for Polar requirements
    public static double getPolarReq(String AminoAcid) {
        switch (AminoAcid) {
            case "Ala":
                return 7.0;
            case "Asp":
                return 13.0;
            case "Cys":
                return 4.8;
            case "Gln":
                return 8.6;
            case "His":
                return 8.4;
            case "Leu":
                return 4.9;
            case "Met":
                return 5.3;
            case "Pro":
                return 6.6;
            case "Thr":
                return 6.6;
            case "Tyr":
                return 5.4;
            case "Arg":
                return 9.1;
            case "Asn":
                return 10.0;
            case "Glu":
                return 12.5;
            case "Gly":
                return 7.9;
            case "Ile":
                return 4.9;
            case "Lys":
                return 10.1;
            case "Phe":
                return 5.0;
            case "Ser":
                return 7.5;
            case "Trp":
                return 5.2;
            case "Val":
                return 5.6;
        }
        throw new InvalidParameterException("Invalid Amino Acid: " + AminoAcid);
    }

    public static int getNucleotideNumber(char chr) throws IllegalArgumentException {
        switch (chr) {
            case 'T':
                return 0;
            case 'C':
                return 1;
            case 'A':
                return 2;
            case 'G':
                return 3;
            default:
                throw new IllegalArgumentException("only CTAG are valid Nucleotid Characters");
        }
    }


    public static boolean isStopCodon(int a, int b, int c) {
        //T
        if (a == 0) {
            //A
            if (b == 2) {
                //A
                if (c == 2) {
                    return true;
                    //G
                } else if (c == 3) {
                    return false;
                }
                //G
            } else if (b == 3) {
                //A
                if (c == 2) {
                    return true;
                }
            }
        }
        return false;
    }

    public static boolean isStopCodon(String codon) {
        if (codon.equalsIgnoreCase("TAA") || codon.equalsIgnoreCase("TAG") || codon.equalsIgnoreCase("TGA")){
            return true;
        }
        return false;
    }

}
