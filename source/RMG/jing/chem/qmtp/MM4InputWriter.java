package jing.chem.qmtp;

import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Map;

import jing.chem.QMTP;
import jing.chem.molFile;
import jing.rxnSys.Logger;

public class MM4InputWriter extends QMInputWriter implements QMInputWritable {

	public String InChIAug;
	
	private class MM4KEYS{
		public static final String INPUT = "Input";
	}
	
	public MM4InputWriter(
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
		//Step 1: write the script for MM4 batch operation
		//	Example script file:
		//	#! /bin/csh
		//	cp testEthylene.mm4 CPD.MM4
		//	cp $MM4_DATDIR/BLANK.DAT PARA.MM4
		//	cp $MM4_DATDIR/CONST.MM4 .
		//	$MM4_EXEDIR/mm4 <<%
		//	1
		//	2
		//	0
		//	%
		//	mv TAPE4.MM4 testEthyleneBatch.out
		//	mv TAPE9.MM4 testEthyleneBatch.opt
		//	exit
		String inpKeyStr = "";
		try{
			//create batch file with executable permissions: cf. http://java.sun.com/docs/books/tutorial/essential/io/fileAttr.html#posix
			File inpKey = new File(dir+"/"+name+".com");
			inpKeyStr="#! /bin/csh\n";
			inpKeyStr+="cp "+name+".mm4 CPD.MM4\n";
			inpKeyStr+="cp $MM4_DATDIR/BLANK.DAT PARA.MM4\n";
			inpKeyStr+="cp $MM4_DATDIR/CONST.MM4 .\n";
			inpKeyStr+="$MM4_EXEDIR/mm4 <<%\n";
			inpKeyStr+="1\n";//read from first line of .mm4 file
			if (!QMTP.useCanTherm){
				if(attemptNumber%scriptAttempts==1) inpKeyStr+="2\n"; //Block-Diagonal Method then Full-Matrix Method
				else if(attemptNumber%scriptAttempts==0) inpKeyStr+="3\n"; //Full-Matrix Method only
				else throw new Exception();//this point should not be reached
			}
			else{//CanTherm case: write the FORCE.MAT file
				if(attemptNumber%scriptAttempts==1) inpKeyStr+="4\n"; //Block-Diagonal Method then Full-Matrix Method
				else if(attemptNumber%scriptAttempts==0) inpKeyStr+="5\n"; //Full-Matrix Method only
				else throw new Exception();//this point should not be reached
				inpKeyStr+="\n";//<RETURN> for temperature
				inpKeyStr+="4\n";//unofficial option 4 for vibrational eigenvector printout to generate Cartesian force constant matrix in FORCE.MAT file
				inpKeyStr+="0\n";//no vibrational amplitude printout
			}
			inpKeyStr+="0\n";//terminate the job
			inpKeyStr+="%\n";
			inpKeyStr+="mv TAPE4.MM4 "+name+".mm4out\n";
			inpKeyStr+="mv TAPE9.MM4 "+name+".mm4opt\n";
			if(QMTP.useCanTherm){
				inpKeyStr+="mv FORCE.MAT "+name+".fmat\n";
			}
			inpKeyStr+="exit\n";
			FileWriter fw = new FileWriter(inpKey);
			fw.write(inpKeyStr);
			fw.close();
		}
		catch(Exception e){
			String err = "Error in writing MM4 script file\n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}
		keywords = new HashMap<String, String>();
		keywords.put(MM4KEYS.INPUT, inpKeyStr);
		return keywords;

	}

	@Override
	public File createInputFile() {
		//Step 2: call the MoleCoor process to create the MM4 input file from the mole file
		try{
			File runningdir=new File(dir);
			//this will only be run on Linux so we don't have to worry about Linux vs. Windows issues
			String command = "python "+System.getenv("RMG")+"/scripts/MM4InputFileMaker.py ";
			//first argument: input file path; for the first attempts, we will use UFF refined coordinates; if that doesn't work, (e.g. case of cyclopropene, InChI=1/C3H4/c1-2-3-1/h1-2H,3H2 OOXWYYGXTJLWHA-UHFFFAOYAJ) we will try crude coordinates (.cmol suffix)
			if(attemptNumber<=scriptAttempts){
				command=command.concat(molfile.getPath() + " ");
			}
			else{
				command=command.concat(molfile.getCrudePath() + " ");
			}
			//second argument: output path
			String inpfilepath=dir+"/"+name+".mm4";
			command=command.concat(inpfilepath+ " ");
			//third argument: molecule name (the augmented InChI)
			command=command.concat(InChIAug+ " ");
			//fourth argument: PYTHONPATH
			command=command.concat(System.getenv("RMG")+"/source/MoleCoor");//this will pass $RMG/source/MoleCoor to the script (in order to get the appropriate path for importing
			Process molecoorProc = Runtime.getRuntime().exec(command, null, runningdir);

		}
		catch(Exception e){
			String err = "Error in running MoleCoor MOL to .MM4 process \n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}
		return new File(name + ".mm4");
	}

}
