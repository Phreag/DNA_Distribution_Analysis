package keser_master.MultiParam;
import keser_master.MultiParam.code.IGreedyCode;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.stream.Collectors;


public abstract class CodeFinderAbstract {
	IGreedyCode naturalCode;

	List<IGreedyCode> currentCodeSetWithValues;
	List<String[]> codesWithoutValues;

	IGreedyCode[] ValueBuffer;
	int nextValue;
	int Threads=5;

	//runs the code finder algorithm.
	public void RunCodeFinder(int maxIterationCount){
		calculateNaturalCodeValues();
		loadCodesFromFile();
		for (int i=0;i<maxIterationCount;i++){
			calculateValues();
			printStatistics(i+1);
			if(i==maxIterationCount-1)break;
			filterCodes(1000);
			mutateCodes();
			clearDuplicates();
			if(i % 5 == 0 && i > 0){
                writeResultsToFile();
            }
		}
		writeResultsToFile();
        System.out.println("finished");
	}

	abstract void calculateNaturalCodeValues();

	private void loadCodesFromFile(){
        codesWithoutValues = new ArrayList<>();
		try {
			BufferedReader br = new BufferedReader(new FileReader("data/codes.txt"));
			String line;
			while ((line = br.readLine()) != null) {
				codesWithoutValues.add(line.split(" ")[1].split("~"));
			}
			br.close();
		} catch (Exception e) {
			System.out.println("FEHLER: Filereader Error");
			e.printStackTrace();
		}
	}



	//Calculates scores for each code in the set.
	private void calculateValues(){
	    nextValue=0;
		currentCodeSetWithValues = new LinkedList<>();
		ValueBuffer = new IGreedyCode[codesWithoutValues.size()];
		//Start the Calculation Threads
		for (int i=0;i<Threads;i++){
			getThreadedCalculator().start();
			try {
				Thread.sleep(1);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}

		//Write Task to save the Results
		//Sleeps for 10ms if it is too fast
        for (int i=0;i<codesWithoutValues.size();i++){
            if(ValueBuffer[i] != null){
                currentCodeSetWithValues.add(ValueBuffer[i]);
                if(i%1000 == 0){
                    System.out.print("\r Threads: "+ Threads +", Finished: "+i+"/"+ValueBuffer.length);
                }
            }else{
                try {
                    Thread.sleep(100);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
                i--;
            }
        }
        System.out.println("\r Threads: "+ Threads +", Finished: "+ValueBuffer.length+"/"+ValueBuffer.length);
	}

	//Prints the statistics of the code set on the console
	private void printStatistics(int i){
		double mean=0;
		double min=999999999;
		double max=-999999999;
		int cnt=0;
		for (IGreedyCode code: currentCodeSetWithValues){
			double x=code.getGreedy_score();
			mean=mean+x;
			if (x<min)min=x;
			if (x>max)max=x;
			cnt++;
		}

		mean=mean/(double)cnt;
		System.out.println("Minimum: "+min);
		System.out.println("Maximum: "+max);
		System.out.println("Mean Value: "+mean);
		try {
			FileWriter fw = new FileWriter(new File("data/GreedyResults.log"), true);
			if (i==0)fw.write("Date, Run,  Minimum,  Maximum,  Mean "+"\n");
			fw.write(new SimpleDateFormat("dd.MM.yyyy HH:mm:ss").format(new Date())+", "+i+", ");
			fw.write(min+", "+max+", "+mean+"\n");
			fw.close();
		} catch (IOException e) {
			System.out.println("Filewriter Error");
			e.printStackTrace();
		}
		if(mean<min){
			System.out.println("FEHLER: Mean too small");
			System.exit(0);
		}
	}

	//generates permutations of the n best codes
	private void filterCodes(int filterAmount) {
		//Extract best codes
		System.out.println("Extracting best Codes...");
		currentCodeSetWithValues.sort(Comparator.comparing(IGreedyCode::getGreedy_score));

		if (filterAmount > currentCodeSetWithValues.size()) {
			filterAmount = currentCodeSetWithValues.size();
		}
		codesWithoutValues = currentCodeSetWithValues.subList(0, filterAmount - 1).stream().map(IGreedyCode::getCode).collect(Collectors.toList());
	}

	private void mutateCodes(){
		//Permutation
		System.out.println("Permutating the best codes to find even better ones...");
        List<String[]> mutatedCodes = new ArrayList<>(codesWithoutValues);
		for (String[] Original : codesWithoutValues) {
			for (int pos1 = 0; pos1 < 20; pos1++) {
				for (int pos2 = pos1 + 1; pos2 < 20; pos2++) {
					String[] NewCode = new String[20];
					System.arraycopy(Original, 0, NewCode, 0, 20);
					NewCode[pos1] = Original[pos2];
					NewCode[pos2] = Original[pos1];
					mutatedCodes.add(NewCode);
				}
			}
		}
		//Also include 5-Switch Permutations of the best codes to avoid local minimum
		for (String[] Original : codesWithoutValues) {
			String[] Changed = new String[20];
			System.arraycopy(Original, 0, Changed, 0, 20);
			int x = 0;
			while (x < 5) {
				Random ran = new Random();
				int x1 = ran.nextInt(20);
				int x2 = ran.nextInt(20);
				if (x1 != x2) {
					x++;
					String[] NewCode = new String[20];
					System.arraycopy(Changed, 0, NewCode, 0, 20);
					NewCode[x1] = Changed[x2];
					NewCode[x2] = Changed[x1];
					Changed = NewCode;
				}
			}
			mutatedCodes.add(Changed);
		}
		codesWithoutValues = mutatedCodes;
	}

	//Removes duplicates from the set
	private void clearDuplicates(){
		List<String> RawStrings=new ArrayList<>();
		for (String[] rCode : codesWithoutValues){
			RawStrings.add(String.join("~", rCode));
		}
		HashSet<String> duplicateFree=new LinkedHashSet<>(RawStrings);
		System.out.println("Number of Duplicates cleared:"+(RawStrings.size()-duplicateFree.size()));
		System.out.println("Remaining Codes for next Run:"+duplicateFree.size());
		codesWithoutValues.clear();
		for (String code : duplicateFree){
			codesWithoutValues.add(code.split("~"));
		}
	}
	
	synchronized int getNextValue(){
		return nextValue++;
	}

	abstract Thread getThreadedCalculator();

	private void writeResultsToFile() {
		currentCodeSetWithValues.sort(Comparator.comparing(IGreedyCode::getGreedy_score));
		try {
			FileWriter fw =new FileWriter("data/greedyResultsMultiParam.csv",false);
			fw.write("Natural Code:\n");
			int codeNumber = 0;
			writeValueLine(fw, codeNumber++, naturalCode);
            fw.write("Algorithm result (best codes first):\n");
			for (IGreedyCode code : currentCodeSetWithValues){
				writeValueLine(fw, codeNumber++, code);
			}
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	abstract void writeValueLine(FileWriter writer, int Number, IGreedyCode c);
}
