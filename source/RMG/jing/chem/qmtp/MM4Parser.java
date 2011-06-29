package jing.chem.qmtp;

import jing.chem.ChemGraph;

public class MM4Parser extends QMParser {

	/**
	 * if flag is not specified in constructor, 
	 * it is assumed that steric energy should be retrieved
	 * from the .mm4out file
	 */
	Boolean printStericEnergy = false;
	
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
	public MM4Parser(String name, String dir, ChemGraph p_chemGraph,
			boolean printStericEnergy){
		this(name, dir, p_chemGraph);
		this.printStericEnergy = printStericEnergy;
		/*
		 * option to print stericenergy before molar mass 
		 * (this will only be used in useHindRot cases,
		 *  but it is always read in with this function
		 */
		if (printStericEnergy) command = command.concat(" 1");
	}
}
