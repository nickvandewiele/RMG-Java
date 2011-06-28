package jing.chem.qmtp;

import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Map;

import jing.chem.QMTP;
import jing.chem.molFile;
import jing.rxnSys.Logger;

public class CanThermInputWriter extends QMInputWriter implements QMInputWritable {

	public String InChIAug;
	
	private class CANTHERMKEYS {
		public static final String CANTHERM = "Cantherm";
		public static final String ROTOR = "Rotor";
	};
	public CanThermInputWriter(
		
	}

	@Override
	public File write() {
		Map<String, String> inputKeywords = createKeywords();

		File inputFile = createInputFile();
		
		return inputFile;
	}

	@Override
	public Map<String, String> createKeywords() {
	
	}

	@Override
	public File createInputFile() {
		try{
			File canFile=new File(dir+"/"+name+".can");
			FileWriter fw = new FileWriter(canFile);
			fw.write(keywords.get(CANTHERMKEYS.CANTHERM));
			fw.close();
			if(keywords.get(CANTHERMKEYS.ROTOR) !=null){//write the rotor information
				File rotFile=new File(dir+"/"+name+".rotinfo");
				FileWriter fwr = new FileWriter(rotFile);
				fwr.write(CANTHERMKEYS.ROTOR);
				fwr.close();
			}
		}
		catch(Exception e){
			String err = "Error in writing CanTherm input \n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}
	}

}
