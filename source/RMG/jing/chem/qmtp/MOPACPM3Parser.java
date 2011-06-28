package jing.chem.qmtp;

import jing.chem.ChemGraph;
import jing.chem.qmtp.QMParser;

public class MOPACPM3Parser extends QMParser {

	public MOPACPM3Parser(String name, String dir, ChemGraph p_chemGraph){
		super(name, dir, p_chemGraph);
		
		parsingTool = new CCLibParser();
		
		inputFileExtension = ".out";

		scriptFile = "MopacPM3ParsingScript.py";

		executable = "python ";
		

		command=null;
		if (System.getProperty("os.name").toLowerCase().contains("windows")){//special windows case where paths can have spaces and are allowed to be surrounded by quotes
			command = executable + "\"" + System.getProperty("RMG.workingDirectory") + scripts + scriptFile+"\" ";
			String logfilepath = "\""+dir+"/" + name + inputFileExtension + "\"";
			command=command.concat(logfilepath);
			command=command.concat(" \""+ System.getenv("RMG")+"/source\"");//this will pass $RMG/source to the script (in order to get the appropriate path for importing
		}
		else{//non-Windows case
			command = executable + System.getProperty("RMG.workingDirectory") + scripts + scriptFile+" ";
			String logfilepath = dir + "/" + name + inputFileExtension;
			command = command.concat(logfilepath);
			command = command.concat(" "+ System.getenv("RMG") + "/source");//this will pass $RMG/source to the script (in order to get the appropriate path for importing
		}
	}
}
