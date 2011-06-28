package jing.chem.qmtp;

import jing.chem.ChemGraph;

public class GaussianPM3Parser extends QMParser {

	public GaussianPM3Parser (String name, String dir, ChemGraph p_chemGraph){
		super(name, dir, p_chemGraph);

		parsingTool = new CCLibParser();
		
		inputFileExtension = ".log";

		scriptFile = "GaussianPM3ParsingScript.py";

		executable = "python ";


		
		command = null;
		//special windows case where paths can have spaces and are allowed to be surrounded by quotes
		if (System.getProperty("os.name").toLowerCase().contains("windows")){
			command = executable+"\\";
			command += System.getProperty("RMG.workingDirectory")+scripts;
			command += scriptFile+"\\ ";
			String logfilepath="\""+dir+"/"+name+inputFileExtension+"\\";
			command = command.concat(logfilepath);
			//this will pass $RMG/source to the script (in order to get the appropriate path for importing
			command = command.concat(" \""+ System.getenv("RMG")+"/source\"");
		}
		else{//non-Windows case
			command = executable;
			command += System.getProperty("RMG.workingDirectory")+scripts;
			command += scriptFile+" ";
			String logfilepath=dir+"/"+name+inputFileExtension;
			command=command.concat(logfilepath);
			command=command.concat(" "+ System.getenv("RMG")+"/source");//this will pass $RMG/source to the script (in order to get the appropriate path for importing
		}
	}

}
