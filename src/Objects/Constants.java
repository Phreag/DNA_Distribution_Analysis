package Objects;

import java.security.InvalidParameterException;

public class Constants {

    public static final String[] Bases = {"T", "C", "A", "G"};
    public static final double sigmaPolar = 2.448999795835;
    public static final double sigmaHydro = 2.9113398977103;
    public static boolean polarReqEnabled = true;
    public static boolean hydropathyEnabled = false;
    public static boolean normalizeBySigma = false;
    public static double cutOff = 0;


    public static double getSqareDifference(String oldAminoAcid, String newAminoAcid) {
        double difference = 0;
        if (polarReqEnabled) {
            double polar1 = getPolarReq(oldAminoAcid);
            double polar2 = getPolarReq(newAminoAcid);
            double diff = polar1 - polar2;
            if(normalizeBySigma){
                diff = diff /sigmaPolar;
            }
            difference = Math.pow(diff, 2);
        }
        if (hydropathyEnabled) {
            double hydro1 = getHydropathy(oldAminoAcid);
            double hydro2 = getHydropathy(newAminoAcid);
            double diff = hydro1 - hydro2;
            if(normalizeBySigma){
                diff = diff /sigmaHydro;
            }
            difference += Math.pow(diff, 2);
        }
        if(cutOff>0){
            if(difference > cutOff)return 0;
        }
        return difference;
    }

    //Holds data for Polar requirements
    public static double getPolarReq(String aminoAcid) {
        switch (aminoAcid) {
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
        throw new InvalidParameterException("Invalid Amino Acid: " + aminoAcid);
    }

    //Holds data for hydropathy
    public static double getHydropathy(String AminoAcid) {
        switch (AminoAcid) {
            case "Ala":
                return 1.8;
            case "Asp":
                return -3.5;
            case "Cys":
                return 2.5;
            case "Gln":
                return -3.5;
            case "His":
                return -3.2;
            case "Leu":
                return 3.8;
            case "Met":
                return 1.9;
            case "Pro":
                return -1.6;
            case "Thr":
                return -0.7;
            case "Tyr":
                return -1.3;
            case "Arg":
                return -4.5;
            case "Asn":
                return -3.5;
            case "Glu":
                return -3.5;
            case "Gly":
                return -0.4;
            case "Ile":
                return 4.5;
            case "Lys":
                return -3.9;
            case "Phe":
                return 2.8;
            case "Ser":
                return -0.8;
            case "Trp":
                return -0.9;
            case "Val":
                return 4.2;
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
                    return true;
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
        if (codon.equalsIgnoreCase("TAA") ||
                codon.equalsIgnoreCase("TAG") ||
                codon.equalsIgnoreCase("TGA")) {
            return true;
        }
        return false;
    }

}
