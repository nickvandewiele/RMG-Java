package jing.chem.qmtp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

import jing.chem.QMTP;
import jing.rxnSys.Logger;

/**
 * Verifies whether a QM job (externalized) has succesfully completed by 
 * searching for specific keywords in the output files
 * @author nmvdewie
 *
 */
public class QMVerifier {
	
	String name;
	
	String directory;
	
	String InChIaug;
	
	public QMVerifier(String name, String directory, String InChIaug){
		this.name = name;
		this.directory = directory;
		this.InChIaug = InChIaug;
	}
	
	//returns true if a Gaussian file for the given name and directory (.log suffix) exists and indicates successful completion (same criteria as used after calculation runs); terminates if the InChI doesn't match the InChI in the file or if there is no InChI in the file; returns false otherwise
	public boolean successfulGaussianResultExistsQ(){
		//part of the code is taken from runGaussian code above
		//look in the output file to check for the successful termination of the Gaussian calculation
		//failed jobs will contain the a line beginning with " Error termination" near the end of the file
		File file = new File(directory+"/"+name+".log");
		if(file.exists()){//if the file exists, do further checks; otherwise, we will skip to final statement and return false
			int failureFlag=0;//flag (1 or 0) indicating whether the Gaussian job failed
			int InChIMatch=0;//flag (1 or 0) indicating whether the InChI in the file matches InChIaug; this can only be 1 if InChIFound is also 1;
			int InChIFound=0;//flag (1 or 0) indicating whether an InChI was found in the log file
			int InChIPartialMatch=0;//flag (1 or 0) indicating whether the InChI in the log file is a substring of the InChI in RMG's memory
			int completeFlag=0;
			String logFileInChI="";
			try{
				FileReader in = new FileReader(file);
				BufferedReader reader = new BufferedReader(in);
				String line=reader.readLine();
				while(line!=null){
					if (line.startsWith(" Error termination ")) failureFlag=1;
					else if (line.startsWith(" ******")){//also look for imaginary frequencies
						if (line.contains("imaginary frequencies")) failureFlag=1;
					}
					else if(line.startsWith(" InChI=")){
						logFileInChI = line.trim();
						//continue reading lines until a line of dashes is found (in this way, we can read InChIs that span multiple lines)
						line=reader.readLine();
						while (!line.startsWith(" --------")){
							logFileInChI += line.trim();
							line=reader.readLine();
						}
						InChIFound=1;
						if(logFileInChI.equals(InChIaug)) InChIMatch=1;
						else if(InChIaug.startsWith(logFileInChI)) InChIPartialMatch=1;
					}
					else if (line.startsWith(" Normal termination of Gaussian 03 at")){//check for completion notice at end of file
						if (reader.readLine()==null){//the above line will occur once in the middle of most files (except for monatomic species), but the following line is non-critical (e.g. "Link1:  Proceeding to internal job step number  2.") so reading over this line should not cause a problem
							completeFlag=1;
						}
					}
					line=reader.readLine();
				}
				reader.close();
				in.close();
			}
			catch(Exception e){
				String err = "Error in reading preexisting Gaussian log file \n";
				err += e.toString();
				Logger.logStackTrace(e);
				System.exit(0);
			}
			//if the complete flag is still 0, the process did not complete and is a failure
			if (completeFlag==0) failureFlag=1;
			//if the failure flag is still 0, the process should have been successful
			if (failureFlag==0&&InChIMatch==1){
				Logger.info("Pre-existing successful quantum result for " + name + " ("+InChIaug+") has been found. This log file will be used.");
				return true;
			}
			else if (InChIFound==1 && InChIMatch == 0){//InChIs do not match (most likely due to limited name length mirrored in log file (79 characters), but possibly due to a collision)
				if(InChIPartialMatch == 1){//case where the InChI in memory begins with the InChI in the log file; we will continue and check the input file, printing a warning if there is no match
					File inputFile = new File(directory+"/"+name+".gjf");
					if(inputFile.exists()){//read the Gaussian inputFile
						String inputFileInChI="";
						try{
							FileReader inI = new FileReader(inputFile);
							BufferedReader readerI = new BufferedReader(inI);
							String lineI=readerI.readLine();
							while(lineI!=null){
								if(lineI.startsWith(" InChI=")){
									inputFileInChI = lineI.trim();
								}
								lineI=readerI.readLine();
							}
							readerI.close();
							inI.close();
						}
						catch(Exception e){
							String err = "Error in reading preexisting Gaussian gjf file \n";
							err += e.toString();
							Logger.logStackTrace(e);
							System.exit(0);
						}
						if(inputFileInChI.equals(InChIaug)){
							if(failureFlag==0){
								Logger.info("Pre-existing successful quantum result for " + name + " ("+InChIaug+") has been found. This log file will be used. *Note that input file was read to confirm lack of InChIKey collision (InChI probably more than 79 characters)");
								return true;
							}
							else{//otherwise, failureFlag==1
								Logger.info("Pre-existing quantum result for " + name + " ("+InChIaug+") has been found, but the result was apparently unsuccessful. The file will be overwritten with a new calculation. *Note that input file was read to confirm lack of InChIKey collision (InChI probably more than 79 characters)");
								return false;
							}
						}
						else{
							if(inputFileInChI.equals("")){//InChI was not found in input file
								Logger.info("*****Warning: potential InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " Log file Augmented InChI = "+logFileInChI + " . InChI could not be found in the Gaussian input file. You should manually check that the log file contains the intended species.");
								return true;
							}
							else{//InChI was found but doesn't match
								Logger.critical("Congratulations! You appear to have discovered the first recorded instance of an InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " Gaussian input file Augmented InChI = "+inputFileInChI);
								System.exit(0);
							}
						}
					}
					else{
						Logger.info("*****Warning: potential InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " Log file Augmented InChI = "+logFileInChI + " . Gaussian input file could not be found to check full InChI. You should manually check that the log file contains the intended species.");
						return true;
					}
				}
				else{
					Logger.critical("Congratulations! You appear to have discovered the first recorded instance of an InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " Log file Augmented InChI = "+logFileInChI);
					System.exit(0);
				}
			}
			else if (InChIFound==0){
				Logger.critical("An InChI was not found in file: " +name+".log");
				System.exit(0);
			}
			else if (failureFlag==1){//note these should cover all possible results for this block, and if the file.exists block is entered, it should return from within the block and should not reach the return statement below
				Logger.info("Pre-existing quantum result for " + name + " ("+InChIaug+") has been found, but the result was apparently unsuccessful. The file will be overwritten with a new calculation.");
				return false;
			}
		}
		//we could print a line here for cases where the file doesn't exist, but this would probably be too verbose
		return false;
	}
	//returns true if a successful result exists (either Gaussian or MOPAC)
	//    public boolean [] successfulResultExistsQ(String name, String directory, String InChIaug){
	//        boolean gaussianResult=successfulGaussianResultExistsQ(name, directory, InChIaug);
	//        boolean mopacResult=successfulMOPACResultExistsQ(name, directory, InChIaug);
	//        return (gaussianResult || mopacResult);// returns true if either a successful Gaussian or MOPAC result exists
	//    }

	//returns true if a MOPAC output file for the given name and directory (.out suffix) exists and indicates successful completion (same criteria as used after calculation runs); terminates if the InChI doesn't match the InChI in the file or if there is no InChI in the file; returns false otherwise
	public boolean successfulMopacResultExistsQ(){
		//part of the code is taken from analogous code for Gaussian
		//look in the output file to check for the successful termination of the calculation (assumed to be successful if "description of vibrations appears)
		File file = new File(directory+"/"+name+".out");
		if(file.exists()){//if the file exists, do further checks; otherwise, we will skip to final statement and return false
			int failureFlag=1;//flag (1 or 0) indicating whether the MOPAC job failed
			int failureOverrideFlag=0;//flag (1 or 0) to override success as measured by failureFlag
			int InChIMatch=0;//flag (1 or 0) indicating whether the InChI in the file matches InChIaug; this can only be 1 if InChIFound is also 1;
			int InChIFound=0;//flag (1 or 0) indicating whether an InChI was found in the log file
			int InChIPartialMatch=0;//flag (1 or 0) indicating whether the InChI in the log file is a substring of the InChI in RMG's memory
			String logFileInChI="";
			try{
				FileReader in = new FileReader(file);
				BufferedReader reader = new BufferedReader(in);
				String line=reader.readLine();
				while(line!=null){
					String trimLine= line.trim();
					if (trimLine.equals("DESCRIPTION OF VIBRATIONS")){
						// if(!MopacFileContainsNegativeFreqsQ(name, directory)) failureFlag=0;//check for this line; if it is here, check for negative frequencies
						failureFlag = 0;
					}
					//negative frequencies notice example:
					//         NOTE: SYSTEM IS NOT A GROUND STATE, THEREFORE ZERO POINT
					//         ENERGY IS NOT MEANINGFULL. ZERO POINT ENERGY PRINTED
					//         DOES NOT INCLUDE THE  2 IMAGINARY FREQUENCIES
					else if (trimLine.endsWith("IMAGINARY FREQUENCIES")){
						//  Logger.info("*****Imaginary freqencies found:");
						failureOverrideFlag=1;
					}
					else if (trimLine.equals("EXCESS NUMBER OF OPTIMIZATION CYCLES")){//exceeding max cycles error
						failureOverrideFlag=1;
					}
					else if (trimLine.equals("NOT ENOUGH TIME FOR ANOTHER CYCLE")){//timeout error
						failureOverrideFlag=1;
					}
					else if(line.startsWith(" InChI=")){
						logFileInChI = line.trim();//output files should take up to 240 characters of the name in the input file
						InChIFound=1;
						if(logFileInChI.equals(InChIaug)) InChIMatch=1;
						else if(InChIaug.startsWith(logFileInChI)) InChIPartialMatch=1;
					}
					line=reader.readLine();
				}
				reader.close();
				in.close();
			}
			catch(Exception e){
				String err = "Error in reading preexisting MOPAC output file \n";
				err += e.toString();
				Logger.logStackTrace(e);
				System.exit(0);
			}
			if(failureOverrideFlag==1) failureFlag=1; //job will be considered a failure if there are imaginary frequencies or if job terminates to to excess time/cycles
			//if the failure flag is still 0, the process should have been successful
			if (failureFlag==0&&InChIMatch==1){
				Logger.info("Pre-existing successful MOPAC quantum result for " + name + " ("+InChIaug+") has been found. This log file will be used.");
				return true;
			}
			else if (InChIFound==1 && InChIMatch == 0){//InChIs do not match (most likely due to limited name length mirrored in log file (240 characters), but possibly due to a collision)
				// if(InChIPartialMatch == 1){//case where the InChI in memory begins with the InChI in the log file; we will continue and check the input file, printing a warning if there is no match
				//look in the input file if the InChI doesn't match (apparently, certain characters can be deleted in MOPAC output file for long InChIs)
				File inputFile = new File(directory+"/"+name+".mop");
				if(inputFile.exists()){//read the MOPAC inputFile
					String inputFileInChI="";
					try{
						FileReader inI = new FileReader(inputFile);
						BufferedReader readerI = new BufferedReader(inI);
						String lineI=readerI.readLine();
						while(lineI!=null){
							if(lineI.startsWith("InChI=")){
								inputFileInChI = lineI.trim();
							}
							lineI=readerI.readLine();
						}
						readerI.close();
						inI.close();
					}
					catch(Exception e){
						String err = "Error in reading preexisting MOPAC input file \n";
						err += e.toString();
						Logger.logStackTrace(e);
						System.exit(0);
					}
					if(inputFileInChI.equals(InChIaug)){
						if(failureFlag==0){
							Logger.info("Pre-existing successful MOPAC quantum result for " + name + " ("+InChIaug+") has been found. This log file will be used. *Note that input file was read to confirm lack of InChIKey collision (InChI probably more than 240 characters or characters probably deleted from InChI in .out file)");
							return true;
						}
						else{//otherwise, failureFlag==1
							Logger.info("Pre-existing MOPAC quantum result for " + name + " ("+InChIaug+") has been found, but the result was apparently unsuccessful. The file will be overwritten with a new calculation or Gaussian result (if available) will be used. *Note that input file was read to confirm lack of InChIKey collision (InChI probably more than 240 characters or characters probably deleted from InChI in .out file)");
							return false;
						}
					}
					else{
						if(inputFileInChI.equals("")){//InChI was not found in input file
							Logger.info("*****Warning: potential InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " Log file Augmented InChI = "+logFileInChI + " . InChI could not be found in the MOPAC input file. You should manually check that the output file contains the intended species.");
							return true;
						}
						else{//InChI was found but doesn't match
							Logger.critical("Congratulations! You appear to have discovered the first recorded instance of an InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " MOPAC input file Augmented InChI = " + inputFileInChI + " Log file Augmented InChI = "+logFileInChI);
							System.exit(0);
						}
					}
				}
				else{
					Logger.info("*****Warning: potential InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " Log file Augmented InChI = "+logFileInChI + " . MOPAC input file could not be found to check full InChI. You should manually check that the log file contains the intended species.");
					return true;
				}
				//  }
				//                else{
				//                    Logger.critical("Congratulations! You appear to have discovered the first recorded instance of an InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " MOPAC output file Augmented InChI = "+logFileInChI);
				//                    System.exit(0);
				//                }
			}
			else if (InChIFound==0){
				Logger.critical("An InChI was not found in file: " +name+".out");
				System.exit(0);
			}
			else if (failureFlag==1){//note these should cover all possible results for this block, and if the file.exists block is entered, it should return from within the block and should not reach the return statement below
				Logger.info("Pre-existing MOPAC quantum result for " + name + " ("+InChIaug+") has been found, but the result was apparently unsuccessful. The file will be overwritten with a new calculation or Gaussian result (if available) will be used.");
				return false;
			}
		}
		//we could print a line here for cases where the file doesn't exist, but this would probably be too verbose
		return false;
	}
	//returns true if an MM4 output file for the given name and directory (.mm4out suffix) exists and indicates successful completion (same criteria as used after calculation runs); terminates if the InChI doesn't match the InChI in the file or if there is no InChI in the file; returns false otherwise
	public boolean successfulMM4ResultExistsQ(){
		//part of the code is taken from analogous code for MOPAC (first ~half) and Gaussian (second ~half)
		//look in the output file to check for the successful termination of the calculation (assumed to be successful if "description of vibrations appears)
		int failureFlag=1;//flag (1 or 0) indicating whether the MM4 job failed
		int failureOverrideFlag=0;//flag (1 or 0) to override success as measured by failureFlag
		File file = new File(directory+"/"+name+".mm4out");
		File canFile = new File(directory+"/"+name+".canout");
		int InChIMatch=0;//flag (1 or 0) indicating whether the InChI in the file matches InChIaug; this can only be 1 if InChIFound is also 1;
		int InChIFound=0;//flag (1 or 0) indicating whether an InChI was found in the log file
		int InChIPartialMatch=0;//flag (1 or 0) indicating whether the InChI in the log file is a substring of the InChI in RMG's memory
		if(QMTP.useCanTherm){//if we are using CanTherm, check whether a CanTherm output file exists...if it does, we will continue on, otherwise, we will rerun calculations (including MM4 calculation) from scratch to ensure our atom numbering is consistent; note: if the .canout file exists, we still want to check for the actual MM4 file even if we are using CanTherm and reading CanTherm output because 1. it ensures we have the correct species and don't have an InChI collision 2. it is needed for getting the geometry which is needed for the symmetry number corrections applied to CanTherm output (which doesn't include symmetry number considerations)
			if(!canFile.exists()) return false;
		}
		if(file.exists()){//if the file exists, do further checks; otherwise, we will skip to final statement and return false
			String logFileInChI="";
			try{
				FileReader in = new FileReader(file);
				BufferedReader reader = new BufferedReader(in);
				String line=reader.readLine();
				while(line!=null){
					String trimLine= line.trim();
					if (trimLine.equals("STATISTICAL THERMODYNAMICS ANALYSIS")){
						failureFlag = 0;
					}
					else if (trimLine.endsWith("imaginary frequencies,")){//read the number of imaginary frequencies and make sure it is zero
						String[] split = trimLine.split("\\s+");
						if (Integer.parseInt(split[3])>0){
							Logger.info("*****Imaginary freqencies found:");
							failureOverrideFlag=1;
						}
					}
					else if (trimLine.contains("             0.0     (fir )")){
						if (QMTP.useCanTherm){//zero frequencies are only acceptable when CanTherm is used
							Logger.info("*****Warning: zero freqencies found (values lower than 7.7 cm^-1 are rounded to zero in MM4 output); CanTherm should hopefully correct this:");
						}
						else{
							Logger.info("*****Zero freqencies found:");
							failureOverrideFlag=1;
						}
					}
					else if(trimLine.startsWith("InChI=")){
						logFileInChI = line.trim();//output files should take up to about 60 (?) characters of the name in the input file
						InChIFound=1;
						if(logFileInChI.equals(InChIaug)) InChIMatch=1;
						else if(InChIaug.startsWith(logFileInChI)) InChIPartialMatch=1;
					}
					line=reader.readLine();
				}
				reader.close();
				in.close();
			}
			catch(Exception e){
				String err = "Error in reading preexisting MM4 output file \n";
				err += e.toString();
				Logger.logStackTrace(e);
				System.exit(0);
			}
			if(failureOverrideFlag==1) failureFlag=1; //job will be considered a failure if there are imaginary frequencies or if job terminates to to excess time/cycles
			//if the failure flag is still 0, the process should have been successful
			if (failureFlag==0&&InChIMatch==1){
				Logger.info("Pre-existing successful MM4 result for " + name + " ("+InChIaug+") has been found. This log file will be used.");
				return true;
			}
			else if (InChIFound==1 && InChIMatch == 0){//InChIs do not match (most likely due to limited name length mirrored in log file (79 characters), but possibly due to a collision)
				if(InChIPartialMatch == 1){//case where the InChI in memory begins with the InChI in the log file; we will continue and check the input file, printing a warning if there is no match
					File inputFile = new File(directory+"/"+name+".mm4");
					if(inputFile.exists()){//read the MM4 inputFile
						String inputFileInChI="";
						try{
							FileReader inI = new FileReader(inputFile);
							BufferedReader readerI = new BufferedReader(inI);
							String lineI=readerI.readLine();
							//InChI should be repeated after in the first line of the input file
							inputFileInChI = lineI.trim().substring(80);//extract the string starting with character 81
							readerI.close();
							inI.close();
						}
						catch(Exception e){
							String err = "Error in reading preexisting MM4 .mm4 file \n";
							err += e.toString();
							Logger.logStackTrace(e);
							System.exit(0);
						}
						if(inputFileInChI.equals(InChIaug)){
							if(failureFlag==0){
								Logger.info("Pre-existing successful MM4 result for " + name + " ("+InChIaug+") has been found. This log file will be used. *Note that input file was read to confirm lack of InChIKey collision (InChI probably more than ~60 characters)");
								return true;
							}
							else{//otherwise, failureFlag==1
								Logger.info("Pre-existing MM4 result for " + name + " ("+InChIaug+") has been found, but the result was apparently unsuccessful. The file will be overwritten with a new calculation. *Note that input file was read to confirm lack of InChIKey collision (InChI probably more than ~60 characters)");
								return false;
							}
						}
						else{
							if(inputFileInChI.equals("")){//InChI was not found in input file
								Logger.info("*****Warning: potential InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " Log file Augmented InChI = "+logFileInChI + " . InChI could not be found in the MM4 input file. You should manually check that the log file contains the intended species.");
								return true;
							}
							else{//InChI was found but doesn't match
								Logger.critical("Congratulations! You appear to have discovered the first recorded instance of an InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " MM4 input file Augmented InChI = "+inputFileInChI);
								System.exit(0);
							}
						}
					}
					else{
						Logger.info("*****Warning: potential InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " Log file Augmented InChI = "+logFileInChI + " . MM4 input file could not be found to check full InChI. You should manually check that the log file contains the intended species.");
						return true;
					}
				}
				else{
					Logger.critical("Congratulations! You appear to have discovered the first recorded instance of an InChIKey collision: InChIKey(augmented) = " + name + " RMG Augmented InChI = "+ InChIaug + " Log file Augmented InChI = "+logFileInChI);
					System.exit(0);
				}
			}
			else if (InChIFound==0){
				Logger.critical("An InChI was not found in file: " +name+".mm4out");
				System.exit(0);
			}
			else if (failureFlag==1){//note these should cover all possible results for this block, and if the file.exists block is entered, it should return from within the block and should not reach the return statement below
				Logger.info("Pre-existing MM4 result for " + name + " ("+InChIaug+") has been found, but the result was apparently unsuccessful. The file will be overwritten with a new calculation.");
				return false;
			}
		}
		//we could print a line here for cases where the file doesn't exist, but this would probably be too verbose
		return false;
	}
}
