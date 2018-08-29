package keser_master;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.URL;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map.Entry;

import org.apache.commons.io.FileUtils;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;

public class GenBankConnection {
	
	public List<DNASequence> LoadMixedFileReadingframe() {
		File f=new File("data/HomoSapiens_CodingSequences.fasta");
		System.out.println("Loading file Mixed.fasta to memory ("+(f.length()/1024)+" KB)...");
		try{
			LinkedHashMap<String, DNASequence> a = FastaReaderHelper.readFastaDNASequence(new File("data/HomoSapiens_CodingSequences.fasta"));
			List<DNASequence> Result = new ArrayList<DNASequence>();
			FileWriter fw = new FileWriter(new File("data/SequenceDataUse.log"), true);
			fw.write("\r\n \r\n Results calculated on " + new SimpleDateFormat("dd.MM.yyyy HH:mm:ss").format(new Date())
					+ ":" + "\r\n");
			fw.write("Loaded Sequence filename: data/HomoSapiens_CodingSequences.fasta\r\n");
			fw.write("Usable Sequences:");

			int loadedCount = 0;
			int mod3=0;
			int notMod3=0;
			int notMod3startCodon=0;
			int startCodonSeq=0;
			int startAndEndSeq=0;
			for (Entry<String, DNASequence> entry : a.entrySet()) {

				loadedCount++; // count sequences loaded
				int length = entry.getValue().getSequenceAsString().length();
				System.out.println(
						"Loaded: " + entry.getValue().getOriginalHeader() + " length = " + length + " nucleotides");
				// nur Sequencen die mit Start Codon beginnen und mit Stop Codon
				// enden und den Triplet Raster entsprechen, werden in Result
				// weiterverarrbeitet
				if ((length % 3) == 0) { // ist die Sequenz dur 3 teilbar?
					 mod3++;
					if (entry.getValue().getSequenceAsString().charAt(0) == 'A'
							&& entry.getValue().getSequenceAsString().charAt(1) == 'T'
							&& entry.getValue().getSequenceAsString().charAt(2) == 'G') {
						startCodonSeq++; // beginnt sie mit einem Startcodon?

//						Anpassung Leseraster Stop Codons in der Sequence zulassen
						if (entry.getValue().getSequenceAsString().charAt(length - 3) == 'T') {
							if ((entry.getValue().getSequenceAsString().charAt(length - 2) == 'A')) {
								if ((entry.getValue().getSequenceAsString().charAt(length - 1) == 'A')
										|| (entry.getValue().getSequenceAsString().charAt(length - 1) == 'G')) {
									startAndEndSeq++;
									fw.write(entry.getValue().getOriginalHeader() + "\r\n");
									Result.add(entry.getValue());
								}
							} else if (entry.getValue().getSequenceAsString().charAt(length - 2) == 'G'
									&& entry.getValue().getSequenceAsString().charAt(length - 1) == 'A') {
								startAndEndSeq++;
								fw.write(entry.getValue().getOriginalHeader() + "\r\n");
								Result.add(entry.getValue());

							}
						}

					} 

				}
				if ((length % 3) != 0) {
					notMod3++;
					if (entry.getValue().getSequenceAsString().charAt(0) == 'A'
							&& entry.getValue().getSequenceAsString().charAt(1) == 'T'
							&& entry.getValue().getSequenceAsString().charAt(2) == 'G') {
						notMod3startCodon++;
					}
				}
			}

			fw.write("\r\n Loaded Sequence Count        : " + loadedCount);
			fw.write("\r\n Count not Mod 3 Count : " + notMod3);
			fw.write("\r\n Count  Sequences with not mod 3 but Start : " + notMod3startCodon);
			fw.write("\r\n Count  Mod 3 Count : " + mod3);
			int check = startCodonSeq - startAndEndSeq;
			int ranStart = mod3 - startCodonSeq;
			fw.write("\r\n Count  Sequence (mod3) with random Start : " + ranStart);
			fw.write("\r\n Count  Sequence (mod3) with Start : " + startCodonSeq);
			fw.write("\r\n Count  Sequence (mod3) with Start but without End : " + check);
			fw.write("\r\n Count  Sequence (mod3) with Start/End : " + startAndEndSeq);
			
			fw.close();

			if (Result.size() == 0) {
				return null;
			}

			return Result;

		} catch (Exception e) {
			System.out.println("ERROR: FASTA-File could not be loaded.");
			System.out.println(e.getMessage());
			e.printStackTrace();
			return null;
		}

	}
	
	//Allows loading of files containing multiple sequences
	public List<DNASequence> LoadMixedFile(){
		File f=new File("data/HomoSapiens_CodingSequences.fasta");
		System.out.println("Loading file Mixed.fasta to memory ("+(f.length()/1024)+" KB)...");
		try{
			LinkedHashMap<String, DNASequence> a = FastaReaderHelper.readFastaDNASequence(new File("data/HomoSapiens_CodingSequences.fasta"));
			List<DNASequence> Result=new ArrayList<DNASequence>();
			int invalid=0;
			for (  Entry<String, DNASequence> entry : a.entrySet() ) {
				Result.add (entry.getValue());
			}
			if (Result.size()==0)return null;
			System.out.println("Invalid Sequences: "+invalid);
			return Result;
		} catch (Exception e) {
			System.out.println("ERROR: FASTA-File could not be loaded.");
			System.out.println(e.getMessage());
			e.printStackTrace();
			return null;
		}
		
	}
	//Loads a FASTA file from /data
	//If the File doesnt exist it donwloads the sequence from the GenBank.
	//GeneID can be either an GI (outdated) or accession string
		public DNASequence LoadFastaFile (String GeneID){
			File f=new File("data/"+GeneID+".fasta");
			if (!f.exists()){
				System.out.println("Sequence not found. Downloading from GenBank...");
				if(!DownloadFasta(GeneID)){return null;}
				System.out.println("Repairing File and removing disabled characters...");
				try {
					int headerend=0;
					String File=FileUtils.readFileToString(f);
					for (int i=0;i<File.length();i++){
						if (File.charAt(i)=='\n'||File.charAt(i)=='\r'){
							headerend=i;
							break;
						}
					}
					System.out.println("Headerend = "+headerend);
					String Header=File.substring(0, headerend);
					String Seq= File.substring(headerend);
					Seq=Seq.replaceAll("U", "");
					Seq=Seq.replaceAll("R", "");
					Seq=Seq.replaceAll("Y", "");
					Seq=Seq.replaceAll("K", "");
					Seq=Seq.replaceAll("M", "");
					Seq=Seq.replaceAll("S", "");
					Seq=Seq.replaceAll("W", "");
					Seq=Seq.replaceAll("B", "");
					Seq=Seq.replaceAll("D", "");
					Seq=Seq.replaceAll("H", "");
					Seq=Seq.replaceAll("V", "");
					Seq=Seq.replaceAll("N", "");
					System.out.println("Writing changes to file...");
					File=Header+Seq;
					f.delete();
					f.createNewFile();
					PrintWriter out = new PrintWriter("data/"+GeneID+".fasta");
					out.println(File);
					out.close();
					
					
				} catch (IOException e) {
					e.printStackTrace();
					return null;
				}
			}
			
			System.out.println("Loading file "+GeneID+".fasta to memory ("+(f.length()/1024)+" KB)...");
			try{
				LinkedHashMap<String, DNASequence> a = FastaReaderHelper.readFastaDNASequence(new File("data/"+GeneID+".fasta"));
				for (  Entry<String, DNASequence> entry : a.entrySet() ) {
					System.out.println("Loaded: " + entry.getValue().getOriginalHeader() + " length = " + entry.getValue().getSequenceAsString().length() + " nucleotides" );
					return entry.getValue();
				}
			} catch (Exception e) {
				System.out.println("ERROR: FASTA-File could not be loaded.");
				System.out.println(e.getMessage());
				e.printStackTrace();
				return null;
			}
			return null;
		}
		//Uses GI or accession to download the Sequence from GenBank
		private boolean DownloadFasta(String GeneID){
			URL u;
			try {
				u=new URL("https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id="+GeneID+"&&showgi=1");File f=new File("data/"+GeneID+".fasta");
				FileUtils.copyURLToFile(u, f);
				System.out.println("Download finished.");
			} catch (Exception e) {
				System.out.println("ERROR: Could not complete download!");
				System.out.println(e.getMessage());
				e.printStackTrace();
				return false;
			}
			return true;
		}
		
	}
