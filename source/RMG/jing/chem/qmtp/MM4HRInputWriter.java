package jing.chem.qmtp;

import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Map;

import jing.chem.QMTP;
import jing.chem.molFile;
import jing.rxnSys.Logger;

public class MM4HRInputWriter extends QMInputWriter implements QMInputWritable {

	public String InChIAug;
	
	public MM4HRInputWriter(
			String name, 
			String directory,
			molFile p_molfile, 
			int attemptNumber, 
			int multiplicity, String inChIaug
			) {
		
		super(name, directory, p_molfile, attemptNumber, multiplicity);
		
		scriptAttempts = 2;
		
		maxAttemptNumber = 2 * scriptAttempts;
		
		this.InChIAug = inChIaug;
	}

	@Override
	public File write() {
		Map<String, String> inputKeywords = createKeywords();

		File inputFile = createInputFile();
		
		return inputFile;
	}

	@Override
	public Map<String, String> createKeywords() {
		
		keywords = new HashMap<String, String>();
		keywords.put("Input", inpKeyStr);
		return keywords;

	}

	@Override
	public File createInputFile() {
	
		return new File(name + ".mm4");
	}

}
