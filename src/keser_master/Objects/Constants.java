package keser_master.Objects;

import java.security.InvalidParameterException;

public class Constants {

    public static final String[] Bases = {"T", "C", "A", "G"};
    public static final double sigmaPolar = 2.448999795835;
    public static final double sigmaHydro = 2.9113398977103;

    public static boolean polarReqEnabled = true;
    public static boolean hydropathyEnabled = false;
    public static boolean molVolEnabled = false;
    public static boolean molWeightEnabled = false;
    public static boolean pKaEnabled = false;
    public static boolean pKbEnabled = false;
    public static boolean pIEnabled = false;
    public static boolean normalizeBySigma = false;
    public static double cutOff = 0;


    public static double getSqareDifference(String oldAminoAcid, String newAminoAcid) {
        double difference = 0;
        if (polarReqEnabled) {
            double polar1 = getPolarReq(oldAminoAcid);
            double polar2 = getPolarReq(newAminoAcid);
            double diff = polar1 - polar2;
            if (normalizeBySigma) {
                diff = diff / sigmaPolar;
            }
            difference = Math.pow(diff, 2);
        }
        if (hydropathyEnabled) {
            double hydro1 = getHydropathy(oldAminoAcid);
            double hydro2 = getHydropathy(newAminoAcid);
            double diff = hydro1 - hydro2;
            if (normalizeBySigma) {
                diff = diff / sigmaHydro;
            }
            difference += Math.pow(diff, 2);
        }
        if (molVolEnabled) {
            double value1 = getMolVol(oldAminoAcid);
            double value2 = getMolVol(newAminoAcid);
            double diff = value1 - value2;
            difference += Math.pow(diff, 2);
        }
        if (molWeightEnabled) {
            double value1 = getMolWeight(oldAminoAcid);
            double value2 = getMolWeight(newAminoAcid);
            double diff = value1 - value2;
            difference += Math.pow(diff, 2);
        }
        if (pKaEnabled) {
            double value1 = getpKa(oldAminoAcid);
            double value2 = getpKa(newAminoAcid);
            double diff = value1 - value2;
            difference += Math.pow(diff, 2);
        }
        if (pKbEnabled) {
            double value1 = getpKb(oldAminoAcid);
            double value2 = getpKb(newAminoAcid);
            double diff = value1 - value2;
            difference += Math.pow(diff, 2);
        }
        if (pIEnabled) {
            double value1 = getpI(oldAminoAcid);
            double value2 = getpI(newAminoAcid);
            double diff = value1 - value2;
            difference += Math.pow(diff, 2);
        }
        if (cutOff > 0) {
            if (difference > cutOff) return 0;
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

    //Holds data for molecular Weight
    public static double getMolWeight(String AminoAcid) {
        switch (AminoAcid) {
            case "Ala":
                return 89.09;
            case "Asp":
                return 133.10;
            case "Cys":
                return 121.16;
            case "Gln":
                return 146.15;
            case "His":
                return 155.16;
            case "Leu":
                return 131.17;
            case "Met":
                return 149.21;
            case "Pro":
                return 115.13;
            case "Thr":
                return 119.12;
            case "Tyr":
                return 181.19;
            case "Arg":
                return 174.20;
            case "Asn":
                return 132.12;
            case "Glu":
                return 147.13;
            case "Gly":
                return 75.07;
            case "Ile":
                return 131.17;
            case "Lys":
                return 146.19;
            case "Phe":
                return 165.19;
            case "Ser":
                return 105.09;
            case "Trp":
                return 204.23;
            case "Val":
                return 117.15;
        }
        throw new InvalidParameterException("Invalid Amino Acid: " + AminoAcid);
    }

    //Holds data for molecular Weight
    public static double getpKa(String AminoAcid) {
        switch (AminoAcid) {
            case "Ala":
                return 2.33;
            case "Asp":
                return 1.95;
            case "Cys":
                return 1.91;
            case "Gln":
                return 2.18;
            case "His":
                return 1.70;
            case "Leu":
                return 2.32;
            case "Met":
                return 2.16;
            case "Pro":
                return 1.95;
            case "Thr":
                return 2.20;
            case "Tyr":
                return 2.24;
            case "Arg":
                return 2.03;
            case "Asn":
                return 2.16;
            case "Glu":
                return 2.16;
            case "Gly":
                return 2.34;
            case "Ile":
                return 2.26;
            case "Lys":
                return 2.15;
            case "Phe":
                return 2.18;
            case "Ser":
                return 2.13;
            case "Trp":
                return 2.38;
            case "Val":
                return 2.27;
        }
        throw new InvalidParameterException("Invalid Amino Acid: " + AminoAcid);
    }

    //Holds data for molecular Weight
    public static double getpKb(String AminoAcid) {
        switch (AminoAcid) {
            case "Ala":
                return 9.71;
            case "Asp":
                return 9.66;
            case "Cys":
                return 10.28;
            case "Gln":
                return 9.00;
            case "His":
                return 9.09;
            case "Leu":
                return 9.58;
            case "Met":
                return 9.08;
            case "Pro":
                return 10.47;
            case "Thr":
                return 8.96;
            case "Tyr":
                return 9.04;
            case "Arg":
                return 9.00;
            case "Asn":
                return 8.73;
            case "Glu":
                return 9.58;
            case "Gly":
                return 9.58;
            case "Ile":
                return 9.60;
            case "Lys":
                return 9.16;
            case "Phe":
                return 9.09;
            case "Ser":
                return 9.05;
            case "Trp":
                return 9.34;
            case "Val":
                return 9.52;
        }
        throw new InvalidParameterException("Invalid Amino Acid: " + AminoAcid);
    }

    //Holds data for molecular Weight
    public static double getpI(String AminoAcid) {
        switch (AminoAcid) {
            case "Ala":
                return 6.00;
            case "Asp":
                return 2.77;
            case "Cys":
                return 5.07;
            case "Gln":
                return 5.65;
            case "His":
                return 7.59;
            case "Leu":
                return 5.98;
            case "Met":
                return 5.74;
            case "Pro":
                return 6.30;
            case "Thr":
                return 5.60;
            case "Tyr":
                return 5.66;
            case "Arg":
                return 10.76;
            case "Asn":
                return 5.41;
            case "Glu":
                return 3.22;
            case "Gly":
                return 5.97;
            case "Ile":
                return 6.02;
            case "Lys":
                return 9.74;
            case "Phe":
                return 5.48;
            case "Ser":
                return 5.68;
            case "Trp":
                return 5.89;
            case "Val":
                return 5.96;
        }
        throw new InvalidParameterException("Invalid Amino Acid: " + AminoAcid);
    }

    //Holds data for molecular Weight
    public static double getMolVol(String AminoAcid) {
        switch (AminoAcid) {
            case "Ala":
                return 31;
            case "Asp":
                return 54;
            case "Cys":
                return 55;
            case "Gln":
                return 85;
            case "His":
                return 96;
            case "Leu":
                return 111;
            case "Met":
                return 105;
            case "Pro":
                return 32.5;
            case "Thr":
                return 61;
            case "Tyr":
                return 136;
            case "Arg":
                return 124;
            case "Asn":
                return 56;
            case "Glu":
                return 83;
            case "Gly":
                return 3;
            case "Ile":
                return 111;
            case "Lys":
                return 119;
            case "Phe":
                return 132;
            case "Ser":
                return 32;
            case "Trp":
                return 170;
            case "Val":
                return 84;
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
