package jing.chem.qmtp;

import jing.chem.ChemGraph;

public class MM4Parser extends QMParser {

	
	public MM4Parser(String name, String dir, ChemGraph p_chemGraph){
		super(name, dir, p_chemGraph);
		
		parsingTool = new CCLibParser();
		
		inputFileExtension = ".mm4out";

		scriptFile = "MM4ParsingScript.py";

		executable = "python ";
		
		command = executable 
		+ System.getProperty("RMG.workingDirectory") + scripts + scriptFile;
		
		String logfilepath = dir + "/" + name + inputFileExtension;
		
		command=command.concat(logfilepath);
		
		//this will pass $RMG/source to the script (in order to get the appropriate path for importing
		command = command.concat(" "+ System.getenv("RMG")+"/source");
	}

}
