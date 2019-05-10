package MultiParam;

import Objects.Constants;
import Objects.Constants_Object;
import Objects.GeneCode;
import keser_master.StabilityCalculator;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class CodeComparationMultiParam {
	//Mapping to read Values from Regine
	String[] RegC={"Phe","Leu","Ile","Met","Val","Ser","Pro","Thr","Ala","Tyr","His","Gln","Asn","Lys","Asp","Glu","Cys","Trp","Arg","Gly"};
	int[] DiffCode={16,   -1,   3,    3,     7,    5,   -5,    0,    4,    8,    -8, -8,    -4,    -4,   -1,   -1,   2,    2,    -14,   -4};
	private int CodeCount;
	private MultiParamCode[] ValueBuffer;
	private List<String> CodeBuffer;
	private int nextValue;
	private FileWriter codes;
	private FileWriter values;
	private Random rnd;
	private int Threads = 3;
	private int ThreadsFinished;
	private class ThreadedGMSCalculator extends Thread {
	    public void run() {
	    	GeneCode g=new GeneCode();
	    	Constants_Object constants = new Constants_Object();
	    	StabilityCalculator S=new StabilityCalculator(g, constants);
	    	while (true){
	    		int currentCode=getNextValue();
	    		
	    		//if (currentCode%10000==0)System.out.println("Generating Code "+currentCode);
	    		if(currentCode>=ValueBuffer.length)break;
				String[] rCode=CodeBuffer.get(currentCode).split(" ")[1].split("~");
				MultiParamCode codeData = new MultiParamCode(rCode);
	    		g.changeCode(rCode);
	    		constants.reset();

    			//0: polarReq
    			codeData.setGMS_Polar(S.getGMS());
				//1: hydropathy
				constants.nextMode();
				codeData.setGMS_Hydro(S.getGMS());
				//2: molVol
				constants.nextMode();
				codeData.setGMS_MolVol(S.getGMS());
				//3: molWeight
				constants.nextMode();
				codeData.setGMS_Mr(S.getGMS());
				//4: pKa
				constants.nextMode();
				codeData.setGMS_pKa(S.getGMS());
				//5: pKb
				constants.nextMode();
				codeData.setGMS_pKb(S.getGMS());
				// 6: pI
				constants.nextMode();
				codeData.setGMS_PI(S.getGMS());

				ValueBuffer[currentCode] = codeData;
	    	}
	    	//System.out.println("Threads finished: "+(++ThreadsFinished));
	    }
	}

	public void loadDefaultcodeSet(){
		try{
			Files.copy(Paths.get("data/codes_defaultset.txt"), Paths.get("data/codes.txt"),StandardCopyOption.REPLACE_EXISTING);
		}catch(Exception e){
			System.out.println("Error loading default code set. Terminating.");
			System.exit(0);
		}
		System.out.println("Default Codeset loaded");
	}

	public void loadCodeSet(String name){
		try{
			Files.copy(Paths.get("data/codes/CodeSet_"+ name+".txt"), Paths.get("data/codes.txt"),StandardCopyOption.REPLACE_EXISTING);
		}catch(Exception e){
			loadCodeSet(name);
		}
		System.out.println("Codeset "+ name + " loaded");
	}

	public void generateCodeSet(String name){
        generateCodeSet(name, true);
	}
	private void generateCodeSet( String name, boolean retry){
		if (retry) {
			File f=new File("data/codes/CodeSet_"+ name+".txt");
			if (f.exists()){
				System.out.println("Codeset "+ name + " already exists. Skipping...");
				return;
			}
		}
		CodeCount = 0;
		File f=new File("data/codes.txt");
		f.delete();
		generateCodes();
		try{
		    File folder = new File("data/codes");
		    if(!folder.exists())folder.mkdir();
			Files.copy(Paths.get("data/codes.txt"), Paths.get("data/codes/CodeSet_"+ name+".txt"),StandardCopyOption.REPLACE_EXISTING);
		}catch(Exception e){
			e.printStackTrace();
			if (retry){
				System.out.println("Error generating code set. Retrying.");
				generateCodeSet(name, false);
			}else{
				System.out.println("Error generating code set. Terminating.");
				System.exit(0);
			}

		}
		System.out.println("Codeset "+ name + " generated");
	}

	public MultiParamCode[] calculateValues(){
		CodeCount=0;
		nextValue=0;
		ThreadsFinished=0;
		//Initialize Writer for Value Output
		CodeBuffer=new ArrayList<>();
		try {
			values=new FileWriter("data/codeValues.csv");
		} catch (IOException e) {
			System.out.println("FileWriter Error!");
			return null;
		}
		System.out.println("Buffering Codes...");
		BufferedReader br;
		try {
			br = new BufferedReader(new FileReader("data/codes.txt"));
			String line = null;
	    	while ((line = br.readLine()) != null) {
	    		CodeBuffer.add(line);
	    	}
	    	br.close();
		} catch (Exception e) {
			System.out.println("FEHLER: Filereader Error");
			return null;
		}
		//Intitalize output Buffer Array
		ValueBuffer=new MultiParamCode[CodeBuffer.size()];
		System.out.println("Start calculation with " +Threads+" Threads...");
		//Start the Calculation Threads
		for (int i=0;i<Threads;i++){
			new ThreadedGMSCalculator().start();
			try {
				Thread.sleep(1);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		//Wait for all Threads
		for (int i=0;i<CodeBuffer.size();i++){
			if(ValueBuffer[i]!=null){
				if (ThreadsFinished==Threads)break;
				//if (i%5000==0)System.out.println("Calculated "+i+" Codes with stabilities");
			}else{
				i--;
				try {
					Thread.sleep(100);
				} catch (InterruptedException e) {}

			}
		}
		try {
			values.close();
		} catch (IOException e) {}
		return ValueBuffer;
	}
	private synchronized int getNextValue(){
		return nextValue++;
	}

	
	//Writes random generated Code line to File
	private void writeCodeLine(String[] rCode){
		CodeCount++;{
		try {
			codes.write(CodeCount+" "+rCode[0]+"~"+rCode[1]+"~"+rCode[2]+"~"+rCode[3]+"~"+rCode[4]+"~"+rCode[5]+"~"+rCode[6]+"~"+rCode[7]+"~"+rCode[8]+"~"+rCode[9]+"~"+rCode[10]+"~"+rCode[11]+"~"+rCode[12]+"~"+rCode[13]+"~"+rCode[14]+"~"+rCode[15]+"~"+rCode[16]+"~"+rCode[17]+"~"+rCode[18]+"~"+rCode[19]+'\n');;
		} catch (IOException e) {
			System.out.println("FileWriter Error on "+CodeCount);
		}
		if (CodeCount%100000==0)System.out.println("Generated "+CodeCount+" Codes");
		}
	}
	
	
	//returns a random generated Code
	private String[] getRandomCode(){
		String[] rCode=new String[20];
		for (int i=0;i<20;i++){
			boolean found=false;
			int acid=0;
			while (found==false){
				acid=rnd.nextInt(20);
				if (!(Contains(rCode,GeneCode.naturalCode[acid],i)))found=true;;
			}
			rCode[i]=GeneCode.naturalCode[acid];
		}
		return rCode;
	}
	private boolean Contains(String[]Snip, String Code, int Level){
		for (int i=0;i<Level;i++){
			if (Snip[i].equals(Code))return true;
		}
		return false;
	}
	
	public void generateCodes(){
		rnd=new Random();
		File f=new File("data/codes.txt");
		if (f.exists()){
			System.out.println("codes.txt existiert bereits.");
			return;
		}
		try {
			codes=new FileWriter("data/codes.txt");
			GeneCode g=new GeneCode();
			//1st is Natural code
	    	writeCodeLine(GeneCode.naturalCode);
	    	
	    	for (int i=0;i<100000;i++){
		    	String[] rCode=getRandomCode();
		    	writeCodeLine(rCode);
		    }
	    	codes.close();
		} catch (IOException e) {
			System.out.println("FileWriter Error!");
			return;
		}
	}
}

