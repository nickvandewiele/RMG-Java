package jing.chem.qmtp;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;

import jing.chem.ChemGraph;
import jing.chem.molFile;
import jing.rxnSys.Logger;

public class ThreeDMolFileCreator {

	String name;
	ChemGraph p_chemGraph;

	public ThreeDMolFileCreator(String name, ChemGraph p_chemGraph) {
		this.name = name;
		this.p_chemGraph = p_chemGraph;
	}
	/**
	 * creates a 3D molFile; for monoatomic species, it just returns the 2D molFile
	 * @return
	 */
	public molFile create() {
		//1. create a 2D file
		//use the absolute path for directory, so we can easily reference from other directories in command-line paths
		//can't use RMG.workingDirectory, since this basically holds the RMG environment variable, not the workingDirectory
		String directory = "2Dmolfiles/";
		File dir=new File(directory);
		directory = dir.getAbsolutePath();
		molFile p_2dfile = new molFile(name, directory, p_chemGraph);
		molFile p_3dfile = new molFile();//it seems this must be initialized, so we initialize to empty object
		//2. convert from 2D to 3D using RDKit if the 2D molfile is for a molecule with 2 or more atoms
		int atoms = p_chemGraph.getAtomNumber();
		if(atoms > 1){
			int distGeomAttempts=1;
			if(atoms > 3){//this check prevents the number of attempts from being negative
				distGeomAttempts = 5*(p_chemGraph.getAtomNumber()-3); //number of conformer attempts is just a linear scaling with molecule size, due to time considerations; in practice, it is probably more like 3^(n-3) or something like that
			}
			p_3dfile = embed3D(p_2dfile, distGeomAttempts);
			return p_3dfile;
		}
		else{
			return p_2dfile;
		}
	}
	//embed a molecule in 3D, using RDKit
	public molFile embed3D(molFile twoDmolFile, int numConfAttempts){
		//convert to 3D MOL file using RDKit script
		int flag=0;
		String directory = "3Dmolfiles/";
		File dir=new File(directory);
		directory = dir.getAbsolutePath();//this uses the absolute path for the directory
		String name = twoDmolFile.getName();
		try{   
			File runningdir=new File(directory);
			String command="";
			if (System.getProperty("os.name").toLowerCase().contains("windows")){//special windows case where paths can have spaces and are allowed to be surrounded by quotes
				command = "python \""+System.getProperty("RMG.workingDirectory")+"/scripts/distGeomScriptMolLowestEnergyConf.py\" ";
				String twoDmolpath=twoDmolFile.getPath();
				command=command.concat("\""+twoDmolpath+"\" ");
				command=command.concat("\""+name+".mol\" ");//this is the target file name; use the same name as the twoDmolFile (but it will be in he 3Dmolfiles folder
				command=command.concat("\""+name+".cmol\" ");//this is the target file name for crude coordinates (corresponding to the minimum energy conformation based on UFF refinement); use the same name as the twoDmolFile (but it will be in he 3Dmolfiles folder) and have suffix .cmol
				command=command.concat(numConfAttempts + " ");
				command=command.concat("\"" + System.getenv("RDBASE")+"\"");//pass the $RDBASE environment variable to the script so it can use the approprate directory when importing rdkit
			}    
			else{//non-Windows case
				command = "python "+System.getProperty("RMG.workingDirectory")+"/scripts/distGeomScriptMolLowestEnergyConf.py ";
				String twoDmolpath=twoDmolFile.getPath();
				command=command.concat(""+twoDmolpath+" ");
				command=command.concat(name+".mol ");//this is the target file name; use the same name as the twoDmolFile (but it will be in he 3Dmolfiles folder
				command=command.concat(name+".cmol ");//this is the target file name for crude coordinates (corresponding to the minimum energy conformation based on UFF refinement); use the same name as the twoDmolFile (but it will be in he 3Dmolfiles folder) and have suffix .cmol
				command=command.concat(numConfAttempts + " ");
				command=command.concat(System.getenv("RDBASE"));//pass the $RDBASE environment variable to the script so it can use the approprate directory when importing rdkit
			}
			Process pythonProc = Runtime.getRuntime().exec(command, null, runningdir);
			String killmsg= "Python process for "+twoDmolFile.getName()+" did not complete within 120 seconds, and the process was killed. File was probably not written.";//message to print if the process times out
			//check for errors and display the error if there is one
			InputStream is = pythonProc.getErrorStream();
			InputStreamReader isr = new InputStreamReader(is);
			BufferedReader br = new BufferedReader(isr);
			String line=null;
			while ( (line = br.readLine()) != null) {
				line = line.trim();
				Logger.error(line);
				flag=1;
			}
			//if there was an error, indicate the file and InChI
			if(flag==1){
				Logger.info("RDKit received error (see above) on " + twoDmolFile.getName()+". File was probably not written.");
			}
			int exitValue = pythonProc.waitFor();
			pythonProc.getInputStream().close();
			pythonProc.getOutputStream().close();
			br.close();
			isr.close();
			is.close();
		}
		catch (Exception e) {
			Logger.logStackTrace(e);
			String err = "Error in running RDKit Python process \n";
			err += e.toString();
			Logger.critical(err);
			System.exit(0);
		}



		// gmagoon 6/3/09 comment out InChI checking for now; in any case, the code will need to be updated, as it is copied from my testing code
		//        //check whether the original InChI is reproduced
		//        if(flag==0){
		//            try{
		//                File f=new File("c:/Python25/"+molfilename);
		//                File newFile= new File("c:/Python25/mol3d.mol");
		//                if(newFile.exists()){
		//                    newFile.delete();//apparently renaming will not work unless target file does not exist (at least on Vista)
		//                }
		//                f.renameTo(newFile);
		//                String command = "c:/Users/User1/Documents/InChI-1/cInChI-1.exe c:/Python25/mol3d.mol inchi3d.inchi /AuxNone /DoNotAddH";//DoNotAddH used to prevent adding Hs to radicals (this would be a problem for current RDKit output which doesn't use M RAD notation)
		//                Process inchiProc = Runtime.getRuntime().exec(command);	
		//               // int exitValue = inchiProc.waitFor();
		//                Thread.sleep(200);//****update: can probably eliminate this using buffered reader
		//                inchiProc.destroy();
		//                
		//                //read output file
		//                File outputFile = new File("inchi3d.inchi");
		//                FileReader fr = new FileReader(outputFile);
		//                BufferedReader br = new BufferedReader(fr);
		//        	String line=null;
		//                String inchi3d=null;
		//                while ( (line = br.readLine()) != null) {
		//                        line = line.trim();
		//                        if(line.startsWith("InChI="))
		//                        {
		//                            inchi3d=line;
		//                        }
		//                }
		//                fr.close();
		//                
		//                //return file to original name:
		//                File f2=new File("c:/Python25/mol3d.mol");
		//                File newFile2= new File("c:/Python25/"+molfilename);
		//                if(newFile2.exists()){
		//                    newFile2.delete();
		//                }
		//                f2.renameTo(newFile2);
		//                
		//                //compare inchi3d with input inchi and print a message if they don't match
		//                if(!inchi3d.equals(inchiString)){
		//                    if(inchi3d.startsWith(inchiString)&&inchiString.length()>10){//second condition ensures 1/C does not match 1/CH4; 6 characters for InChI=, 2 characters for 1/, 2 characters for atom layer
		//                        Logger.info("(probably minor) For File: "+ molfilename+" , 3D InChI (" + inchi3d+") begins with, but does not match original InChI ("+inchiString+"). SMILES string: "+ smilesString);
		//                        
		//                    }
		//                    else{
		//                        Logger.info("For File: "+ molfilename+" , 3D InChI (" + inchi3d+") does not match original InChI ("+inchiString+"). SMILES string: "+ smilesString);
		//                    }
		//                }
		//            }
		//            catch (Exception e) {
		//					Logger.logStackTrace(e);
		//                String err = "Error in running InChI process \n";
		//                err += e.toString();
		//                Logger.critical(err);
		//                System.exit(0);
		//            }
		//        }


		//construct molFile pointer to new file (name will be same as 2D mol file
		return new molFile(name, directory);
	}

}
