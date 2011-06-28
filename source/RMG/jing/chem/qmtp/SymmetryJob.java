package jing.chem.qmtp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import jing.rxnSys.Logger;


/**
 *determine the point group using the SYMMETRY program (http://www.cobalt.chem.ucalgary.ca/ps/symmetry/)
 *required input is a line with number of atoms followed by lines for each atom including atom number and x,y,z coordinates
 *finalTol determines how loose the point group criteria are; values are comparable to those specifed in the GaussView point group interface
 *public String determinePointGroupUsingSYMMETRYProgram(String geom, double finalTol){ 
 *
 */
public class SymmetryJob extends QMJob implements ISymmetryJob {

	IQMData data;

	String inputFile = "input.txt";

	int attemptNumber = 1;
	
	int maxAttemptNumber = 4;
	
	boolean pointGroupFound=false;
	
	public SymmetryJob (IQMData data){
		this.data = data;
	}

	@Override
	public PointGroup calculate() {
		String geom = data.getNumberOfAtoms() + "\n";
		for(int i=0; i < data.getNumberOfAtoms(); i++){
			geom += data.getAtomicNumbers().get(i) 
			+ " "+ data.getThreeDCoords().get(i).getX() 
			+ " " +data.getThreeDCoords().get(i).getY()  
			+ " " +data.getThreeDCoords().get(i).getZ()  
			+ "\n";
		}
		writeInputFile(geom);

		String result = "";

		while (attemptNumber<=maxAttemptNumber && !pointGroupFound){
			//call the program and read the result
			result = "";
			String [] lineArray;

			if (System.getProperty("os.name").toLowerCase().contains("windows")){//the Windows case where the precompiled executable seems to need to be called from a batch script
				if(attemptNumber==1) command = "\""+System.getProperty("RMG.workingDirectory")+"/scripts/symmetryDefault2.bat\" "+qmfolder+ "symminput.txt";//12/1/09 gmagoon: switched to use slightly looser criteria of 0.02 rather than 0.01 to handle methylperoxyl radical result from MOPAC
				else if (attemptNumber==2) command = "\""+System.getProperty("RMG.workingDirectory")+"/scripts/symmetryLoose.bat\" " +qmfolder+ "symminput.txt";//looser criteria (0.1 instead of 0.01) to properly identify C2v group in VBURLMBUVWIEMQ-UHFFFAOYAVmult5 (InChI=1/C3H4O2/c1-3(2,4)5/h1-2H2/mult5) MOPAC result; C2 and sigma were identified with default, but it should be C2 and sigma*2
				else if (attemptNumber==3) command = "\""+System.getProperty("RMG.workingDirectory")+"/scripts/symmetryLoose2.bat\" " +qmfolder+ "symminput.txt";//looser criteria to properly identify D2d group in XXHDHKZTASMVSX-UHFFFAOYAM
				else if (attemptNumber==4){
					Logger.error("*****WARNING****: Using last-resort symmetry estimation options; symmetry may be underestimated");
					command = "\""+System.getProperty("RMG.workingDirectory")+"/scripts/symmetryLastResort.bat\" " +qmfolder+ "symminput.txt";//last resort criteria to avoid crashing (this will likely result in identification of C1 point group)
				}
				else{
					Logger.critical("Invalid attemptNumber: "+ attemptNumber);
					System.exit(0);
				}
			}
			else{//in other (non-Windows) cases, where it is compiled from scratch, we should be able to run this directly
				if(attemptNumber==1) command = System.getProperty("RMG.workingDirectory")+"/bin/SYMMETRY.EXE -final 0.02 " +qmfolder+ "symminput.txt";//12/1/09 gmagoon: switched to use slightly looser criteria of 0.02 rather than 0.01 to handle methylperoxyl radical result from MOPAC
				else if (attemptNumber==2) command = System.getProperty("RMG.workingDirectory")+"/bin/SYMMETRY.EXE -final 0.1 " +qmfolder+ "symminput.txt";//looser criteria (0.1 instead of 0.01) to properly identify C2v group in VBURLMBUVWIEMQ-UHFFFAOYAVmult5 (InChI=1/C3H4O2/c1-3(2,4)5/h1-2H2/mult5) MOPAC result; C2 and sigma were identified with default, but it should be C2 and sigma*2
				else if (attemptNumber==3) command = System.getProperty("RMG.workingDirectory")+"/bin/SYMMETRY.EXE -primary 0.2 -final 0.1 " +qmfolder+ "symminput.txt";//looser criteria to identify D2d group in XXHDHKZTASMVSX-UHFFFAOYAM (InChI=1/C12H16/c1-5-9-10(6-2)12(8-4)11(9)7-3/h5-12H,1-4H2)
				else if (attemptNumber==4){
					Logger.warning("*****WARNING****: Using last-resort symmetry estimation options; symmetry may be underestimated");
					command = System.getProperty("RMG.workingDirectory")+"/bin/SYMMETRY.EXE -final 0.0 " +qmfolder+ "symminput.txt";//last resort criteria to avoid crashing (this will likely result in identification of C1 point group)
				}
				else{
					Logger.critical("Invalid attemptNumber: "+ attemptNumber);
					System.exit(0);
				}
			}

			result = run();

			//check for a recognized point group
			if(PointGroupDictionary.getInstance().contains(result)){
				pointGroupFound=true;
			}
			else{
				if(attemptNumber < maxAttemptNumber) Logger.info("Attempt number "+attemptNumber+" did not identify a recognized point group (" +result+"). Will retry with looser point group criteria.");
				else{
					Logger.critical("Final attempt number "+attemptNumber+" did not identify a recognized point group (" +result+"). Exiting.");
					System.exit(0);
				}
				attemptNumber++;
			}
			
		}
		return new PointGroup(result);

	}
	public String run() {
		Process job = null;
		try{

			job = Runtime.getRuntime().exec(command);
		}
		catch(Exception e){
			String err = "Error in running point group calculation process using SYMMETRY \n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}
		return check(job);
	}


	public String check(Process job) {
		String result = "";
		try{
			//check for errors and display the error if there is one
			InputStream is = job.getInputStream();
			InputStreamReader isr = new InputStreamReader(is);
			BufferedReader br = new BufferedReader(isr);
			String line=null;
			String [] lineArray;
			while ( (line = br.readLine()) != null) {
				if(line.startsWith("It seems to be the ")){//last line, ("It seems to be the [x] point group") indicates point group
					lineArray = line.split(" ");//split the line around spaces
					result = lineArray[5];//point group string should be the 6th word
				}
			}
			int exitValue = job.waitFor();
			job.getErrorStream().close();
			job.getOutputStream().close();
			br.close();
			isr.close();
			is.close();
		}
		catch(Exception e){
			String err = "Error in veryfying point group calculation process using SYMMETRY \n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}
		Logger.info("Point group: "+ result);//print result, at least for debugging purposes
		return result;
	}
	public File writeInputFile(String geom){
		//write the input file
		File path = null;
		try {
			path=new File(qmfolder+inputFile);//SYMMETRY program directory
			FileWriter fw = new FileWriter(path);
			fw.write(geom);
			fw.close();
		} catch (IOException e) {
			String err = "Error writing input file for point group calculation";
			err += e.toString();
			Logger.critical(err);
			System.exit(0);
		}
		return path;
	}
}
