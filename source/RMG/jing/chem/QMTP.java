////////////////////////////////////////////////////////////////////////////////
//
//	RMG - Reaction Mechanism Generator
//
//	Copyright (c) 2002-2011 Prof. William H. Green (whgreen@mit.edu) and the
//	RMG Team (rmg_dev@mit.edu)
//
//	Permission is hereby granted, free of charge, to any person obtaining a
//	copy of this software and associated documentation files (the "Software"),
//	to deal in the Software without restriction, including without limitation
//	the rights to use, copy, modify, merge, publish, distribute, sublicense,
//	and/or sell copies of the Software, and to permit persons to whom the
//	Software is furnished to do so, subject to the following conditions:
//
//	The above copyright notice and this permission notice shall be included in
//	all copies or substantial portions of the Software.
//
//	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//	DEALINGS IN THE SOFTWARE.
//
////////////////////////////////////////////////////////////////////////////////



package jing.chem;



import java.util.*;

import jing.chem.qmtp.CanThermJob;
import jing.chem.qmtp.GaussianJob;
import jing.chem.qmtp.GaussianPM3InputWriter;
import jing.chem.qmtp.GaussianPM3Parser;
import jing.chem.qmtp.IParsingTool;
import jing.chem.qmtp.IQMData;
import jing.chem.qmtp.MM4HRInputWriter;
import jing.chem.qmtp.MM4HRJob;
import jing.chem.qmtp.MM4InputWriter;
import jing.chem.qmtp.MM4Job;
import jing.chem.qmtp.MM4Parser;
import jing.chem.qmtp.MOPACJob;
import jing.chem.qmtp.MOPACPM3InputWriter;
import jing.chem.qmtp.MOPACPM3Parser;
import jing.chem.qmtp.QMConstants;
import jing.chem.qmtp.QMInputWritable;
import jing.chem.qmtp.QMInputWriter;
import jing.chem.qmtp.QMJob;
import jing.chem.qmtp.QMJobRunnable;
import jing.chem.qmtp.QMParsable;
import jing.chem.qmtp.QMParser;
import jing.chemUtil.*;
import jing.param.*;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;

import quicktime.qd3d.math.Point3D;
import jing.rxnSys.Logger;

//quantum mechanics thermo property estimator; analog of GATP
public class QMTP implements GeneralGAPP {

	private static QMTP INSTANCE = new QMTP();		//## attribute INSTANCE
	protected static PrimaryThermoLibrary primaryLibrary;//Note: may be able to separate this out into GeneralGAPP, as this is common to both GATP and QMTP
	public static String qmfolder= "QMfiles/";
	//   protected static HashMap library;		//as above, may be able to move this and associated functions to GeneralGAPP (and possibly change from "x implements y" to "x extends y"), as it is common to both GATP and QMTP
	protected ThermoGAGroupLibrary thermoLibrary; //needed for HBI
	public static String qmprogram= "both";//the qmprogram can be "mopac", "gaussian03", "both" (MOPAC and Gaussian), or "mm4"/"mm4hr"
	public static boolean usePolar = false; //use polar keyword in MOPAC
	public static boolean useCanTherm = true; //whether to use CanTherm in MM4 cases for interpreting output via force-constant matrix; this will hopefully avoid zero frequency issues
	public static boolean useHindRot = false;//whether to use HinderedRotor scans with MM4 (requires useCanTherm=true)
	public static double deltaTheta=5.0;//degree increment for rotor scans when using useHindRot
	// Constructors

	//## operation QMTP()
	private QMTP() {
		// initializeLibrary(); //gmagoon 72509: commented out in GATP, so I am mirroring the change here; other library functions below also commented out
		initializePrimaryThermoLibrary();
	}
	//## operation generateThermoData(ChemGraph)
	public ThermoData generateThermoData(ChemGraph p_chemGraph) {
		//#[ operation generateThermoData(ChemGraph)
		//first, check for thermo data in the primary thermo library and library (?); if it is there, use it
		ThermoData result = primaryLibrary.getThermoData(p_chemGraph.getGraph());
		//Logger.info(result);
		if (result != null) {
			p_chemGraph.fromprimarythermolibrary = true;
			return result;
		}


		//        result = getFromLibrary(p_chemGraph.getChemicalFormula());//gmagoon 72509: commented out in GATP, so I am mirroring the change here
		//        if (result != null) return result;

		result=new ThermoData();

		int maxRadNumForQM = Global.maxRadNumForQM;
		if (p_chemGraph.getRadicalNumber() > maxRadNumForQM)//use HBI if the molecule has more radicals than maxRadNumForQM; this is helpful because ; also MM4 (and MM3) look like they may have issues with radicals
		{//this code is based closely off of GATP saturation (in getGAGroup()), but there are some modifications, particularly for symmetry correction
			//find the initial symmetry number
			int sigmaRadical = p_chemGraph.getSymmetryNumber();

			Graph g = p_chemGraph.getGraph();
			HashMap oldCentralNode = (HashMap)(p_chemGraph.getCentralNode()).clone();
			// saturate radical site
			int max_radNum_molecule = ChemGraph.getMAX_RADICAL_NUM();
			int max_radNum_atom = Math.min(8,max_radNum_molecule);
			int [] idArray = new int[max_radNum_molecule];
			Atom []  atomArray = new Atom[max_radNum_molecule];
			Node [][] newnode = new Node[max_radNum_molecule][max_radNum_atom];

			int radicalSite = 0;
			Iterator iter = p_chemGraph.getNodeList();
			FreeElectron satuated = FreeElectron.make("0");
			while (iter.hasNext()) {
				Node node = (Node)iter.next();
				Atom atom = (Atom)node.getElement();
				if (atom.isRadical()) {
					radicalSite ++;
					// save the old radical atom
					idArray[radicalSite-1] = node.getID().intValue();
					atomArray[radicalSite-1] = atom;
					// new a satuated atom and replace the old one
					Atom newAtom = new Atom(atom.getChemElement(),satuated);
					node.setElement(newAtom);
					node.updateFeElement();
				}
			}

			// add H to saturate chem graph
			Atom H = Atom.make(ChemElement.make("H"),satuated);
			Bond S = Bond.make("S");
			for (int i=0;i<radicalSite;i++) {
				Node node = p_chemGraph.getNodeAt(idArray[i]);
				Atom atom = atomArray[i];
				int HNum = atom.getRadicalNumber();
				for (int j=0;j<HNum;j++) {
					newnode[i][j] = g.addNode(H);
					g.addArcBetween(node,S,newnode[i][j]);
				}
				node.updateFgElement();
			}

			//find the saturated symmetry number
			int sigmaSaturated = p_chemGraph.getSymmetryNumber();

			//         result = generateThermoData(g);//I'm not sure what GATP does, but this recursive calling will use HBIs on saturated species if it exists in PrimaryThermoLibrary
			//check the primary thermo library for the saturated graph
			result = primaryLibrary.getThermoData(p_chemGraph.getGraph());
			//Logger.info(result);
			if (result != null) {
				p_chemGraph.fromprimarythermolibrary = true;
			}
			else{
				result=generateQMThermoData(p_chemGraph);
			}

			// find the BDE for all radical groups
			if(thermoLibrary == null) initGAGroupLibrary();
			for (int i=0; i<radicalSite; i++) {
				int id = idArray[i];
				Node node = g.getNodeAt(id);
				Atom old = (Atom)node.getElement();
				node.setElement(atomArray[i]);
				node.updateFeElement();

				// get rid of the extra H at ith site
				int HNum = atomArray[i].getRadicalNumber();
				for (int j=0;j<HNum;j++) {
					g.removeNode(newnode[i][j]);
				}
				node.updateFgElement();

				p_chemGraph.resetThermoSite(node);
				ThermoGAValue thisGAValue = thermoLibrary.findRadicalGroup(p_chemGraph);
				if (thisGAValue == null) {
					Logger.error("Radical group not found: " + node.getID());
				}
				else {
					//Logger.info(node.getID() + " radical correction: " + thisGAValue.getName() + "  "+thisGAValue.toString());
					result.plus(thisGAValue);
				}

				//recover the saturated site for next radical site calculation
				node.setElement(old);
				node.updateFeElement();
				for (int j=0;j<HNum;j++) {
					newnode[i][j] = g.addNode(H);
					g.addArcBetween(node,S,newnode[i][j]);
				}
				node.updateFgElement();

			}

			// recover the chem graph structure
			// recover the radical
			for (int i=0; i<radicalSite; i++) {
				int id = idArray[i];
				Node node = g.getNodeAt(id);
				node.setElement(atomArray[i]);
				node.updateFeElement();
				int HNum = atomArray[i].getRadicalNumber();
				//get rid of extra H
				for (int j=0;j<HNum;j++) {
					g.removeNode(newnode[i][j]);
				}
				node.updateFgElement();
			}

			// subtract the enthalphy of H from the result
			int rad_number = p_chemGraph.getRadicalNumber();
			ThermoGAValue enthalpy_H = new ThermoGAValue(QMConstants.ENTHALPY_HYDROGEN * rad_number, 0,0,0,0,0,0,0,0,0,0,0,null);
			result.minus(enthalpy_H);


			//correct the symmetry number based on the relative radical and saturated symmetry number; this should hopefully sidestep potential complications based on the fact that certain symmetry effects could be included in HBI value itself, and the fact that the symmetry number correction for saturated molecule has already been implemented, and it is likely to be different than symmetry number considered here, since the correction for the saturated molecule will have been external symmetry number, whereas RMG's ChemGraph symmetry number estimator includes both internal and external symmetry contributions; even so, I don't know if this will handle a change from chiral to achiral (or vice versa) properly            
			ThermoGAValue symmetryNumberCorrection = new ThermoGAValue(0,-1*GasConstant.getCalMolK()*Math.log((double)(sigmaRadical)/(double)(sigmaSaturated)),0,0,0,0,0,0,0,0,0,0,null);
			result.plus(symmetryNumberCorrection);

			p_chemGraph.setCentralNode(oldCentralNode);  

			//display corrected thermo to user
			String [] InChInames = getQMFileName(p_chemGraph);//determine the filename (InChIKey) and InChI with appended info for triplets, etc.
			String name = InChInames[0];
			String InChIaug = InChInames[1];
			Logger.info("HBI-based thermo for " + name + "("+InChIaug+"): "+ result.toString());//print result, at least for debugging purposes
		}
		else{
			result = generateQMThermoData(p_chemGraph);
		}

		return result;
		//#]
	}

	public ThermoData generateQMThermoData(ChemGraph p_chemGraph){
		//if there is no data in the libraries, calculate the result based on QM or MM calculations; the below steps will be generalized later to allow for other quantum mechanics packages, etc.
		String qmProgram = qmprogram;
		String qmMethod = "";
		if(qmProgram.equals("mm4")||qmProgram.equals("mm4hr")){
			qmMethod = "mm4";
			if(qmProgram.equals("mm4hr")) useHindRot=true;
		}
		else{
			qmMethod="pm3"; //may eventually want to pass this to various functions to choose which "sub-function" to call
		}

		ThermoData result = new ThermoData();
		double[] dihedralMinima = null;

		String [] InChInames = getQMFileName(p_chemGraph);//determine the filename (InChIKey) and InChI with appended info for triplets, etc.
		String name = InChInames[0];
		String InChIaug = InChInames[1];
		String directory = qmfolder;
		File dir=new File(directory);
		directory = dir.getAbsolutePath();//this and previous three lines get the absolute path for the directory
		if(qmMethod.equals("pm3")){
			//first, check to see if the result already exists and the job terminated successfully
			boolean gaussianResultExists = successfulGaussianResultExistsQ(name,directory,InChIaug);
			boolean mopacResultExists = successfulMopacResultExistsQ(name,directory,InChIaug);
			if(!gaussianResultExists && !mopacResultExists){//if a successful result doesn't exist from previous run (or from this run), run the calculation; if a successful result exists, we will skip directly to parsing the file
				//steps 1 and 2: create 2D and 3D mole files
				molFile p_3dfile = create3Dmolfile(name, p_chemGraph);
				//3. create the Gaussian or MOPAC input file
				directory = qmfolder;
				dir=new File(directory);
				directory = dir.getAbsolutePath();//this and previous three lines get the absolute path for the directory
				int attemptNumber=1;//counter for attempts using different keywords
				int successFlag=0;//flag for success of Gaussian run; 0 means it failed, 1 means it succeeded
				int maxAttemptNumber=1;
				int multiplicity = p_chemGraph.getRadicalNumber()+1; //multiplicity = radical number + 1
				while(successFlag==0 && attemptNumber <= maxAttemptNumber){
					//IF block to check which program to use
					if (qmProgram.equals("gaussian03")){
						if(p_chemGraph.getAtomNumber() > 1){
							maxAttemptNumber = createGaussianPM3Input(name, directory, p_3dfile, attemptNumber, InChIaug, multiplicity);
						}
						else{
							maxAttemptNumber = createGaussianPM3Input(name, directory, p_3dfile, -1, InChIaug, multiplicity);//use -1 for attemptNumber for monoatomic case
						}
						//4. run Gaussian
						successFlag = runGaussian(name, directory);
					}
					else if (qmProgram.equals("mopac") || qmProgram.equals("both")){
						maxAttemptNumber = createMopacPM3Input(name, directory, p_3dfile, attemptNumber, InChIaug, multiplicity);
						successFlag = runMOPAC(name, directory);
					}
					else{
						Logger.critical("Unsupported quantum chemistry program");
						System.exit(0);
					}
					//new IF block to check success
					if(successFlag==1){
						Logger.info("Attempt #"+attemptNumber + " on species " + name + " ("+InChIaug+") succeeded.");
					}
					else if(successFlag==0){
						if(attemptNumber==maxAttemptNumber){//if this is the last possible attempt, and the calculation fails, exit with an error message
							if(qmProgram.equals("both")){ //if we are running with "both" option and all keywords fail, try with Gaussian
								qmProgram = "gaussian03";
								Logger.info("*****Final MOPAC attempt (#" + maxAttemptNumber + ") on species " + name + " ("+InChIaug+") failed. Trying to use Gaussian.");
								attemptNumber=0;//this needs to be 0 so that when we increment attemptNumber below, it becomes 1 when returning to the beginning of the for loop
								maxAttemptNumber=1;
							}
							else{
								Logger.info("*****Final attempt (#" + maxAttemptNumber + ") on species " + name + " ("+InChIaug+") failed.");
								Logger.critical(p_chemGraph.toString());
								System.exit(0);
								//	ThermoData temp = new ThermoData(1000,0,0,0,0,0,0,0,0,0,0,0,"failed calculation");
								//	temp.setSource("***failed calculation***");
								//	return temp;
							}
						}
						Logger.info("*****Attempt #"+attemptNumber + " on species " + name + " ("+InChIaug+") failed. Will attempt a new keyword.");
						attemptNumber++;//try again with new keyword
					}
				}


			}
			//5. parse QM output and record as thermo data (function includes symmetry/point group calcs, etc.); if both Gaussian and MOPAC results exist, Gaussian result is used
			if (gaussianResultExists || (qmProgram.equals("gaussian03") && !mopacResultExists)){
				result = parseGaussianPM3(name, directory, p_chemGraph);
			}
			else if (mopacResultExists || qmProgram.equals("mopac") || qmProgram.equals("both")){
				result = parseMopacPM3(name, directory, p_chemGraph);
			}
			else{
				Logger.critical("Unexpected situation in QMTP thermo estimation");
				System.exit(0);
			}
		}
		else{//mm4 case
			//first, check to see if the result already exists and the job terminated successfully
			boolean mm4ResultExists = successfulMM4ResultExistsQ(name,directory,InChIaug);
			if(!mm4ResultExists){//if a successful result doesn't exist from previous run (or from this run), run the calculation; if a successful result exists, we will skip directly to parsing the file
				//steps 1 and 2: create 2D and 3D mole files
				molFile p_3dfile = create3Dmolfile(name, p_chemGraph);
				//3. create the MM4 input file
				directory = qmfolder;
				dir=new File(directory);
				directory = dir.getAbsolutePath();//this and previous three lines get the absolute path for the directory
				int attemptNumber=1;//counter for attempts using different keywords
				int successFlag=0;//flag for success of MM4 run; 0 means it failed, 1 means it succeeded
				int maxAttemptNumber=1;
				int multiplicity = p_chemGraph.getRadicalNumber()+1; //multiplicity = radical number + 1
				while(successFlag==0 && attemptNumber <= maxAttemptNumber){
					maxAttemptNumber = createMM4Input(name, directory, p_3dfile, attemptNumber, InChIaug, multiplicity);
					//4. run MM4
					successFlag = runMM4(name, directory);
					//new IF block to check success
					if(successFlag==1){
						Logger.info("Attempt #"+attemptNumber + " on species " + name + " ("+InChIaug+") succeeded.");
						//run rotor calculations if necessary
						int rotors = p_chemGraph.getInternalRotor();
						if(useHindRot && rotors > 0){
							//we should re-run scans even if pre-existing scans exist because atom numbering may change from case to case; a better solution would be to check for stored CanTherm output and use that if available
							Logger.info("Running rotor scans on "+name+"...");
							dihedralMinima = createMM4RotorInput(name, directory, p_chemGraph, rotors);//we don't worry about checking InChI here; if there is an InChI mismatch it should be caught
							runMM4Rotor(name, directory, rotors);
						}
					}
					else if(successFlag==0){
						if(attemptNumber==maxAttemptNumber){//if this is the last possible attempt, and the calculation fails, exit with an error message
							Logger.info("*****Final attempt (#" + maxAttemptNumber + ") on species " + name + " ("+InChIaug+") failed.");
							Logger.critical(p_chemGraph.toString());
							System.exit(0);
							//	ThermoData temp = new ThermoData(1000,0,0,0,0,0,0,0,0,0,0,0,"failed calculation");
							//	temp.setSource("***failed calculation***");
							//	return temp;
						}
						Logger.info("*****Attempt #"+attemptNumber + " on species " + name + " ("+InChIaug+") failed. Will attempt a new keyword.");
						attemptNumber++;//try again with new keyword
					}
				}
				if(useCanTherm){
					performCanThermCalcs(name, directory, p_chemGraph, dihedralMinima, false);
					if (p_chemGraph.getInternalRotor()>0) performCanThermCalcs(name, directory, p_chemGraph, dihedralMinima, true);//calculate RRHO case for comparison
				}

			}
			//5. parse MM4 output and record as thermo data (function includes symmetry/point group calcs, etc.)
			if(!useCanTherm) result = parseMM4(name, directory, p_chemGraph);
			else{
				//if (qmdata==null) qmdata = getQMDataWithCClib(name, directory, p_chemGraph, true);//get qmdata if it is null (i.e. if a pre-existing successful result exists and it wasn't read in above)
				result = parseCanThermFile(name, directory, p_chemGraph);
				if (p_chemGraph.getInternalRotor()>0) parseCanThermFile(name+"_RRHO", directory, p_chemGraph);//print the RRHO result for comparison
			}
		}

		return result;
	}

	protected static QMTP getINSTANCE() {
		return INSTANCE;
	}


	public void initializePrimaryThermoLibrary(){//svp

		primaryLibrary = PrimaryThermoLibrary.getINSTANCE();

	}

	//creates a 3D molFile; for monoatomic species, it just returns the 2D molFile
	public molFile create3Dmolfile(String name, ChemGraph p_chemGraph){
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

	/**
	 *creates Gaussian PM3 input file in directory with filename name.gjf by using OpenBabel
	 * to convert p_molfile
	 *@param attemptNumber determines which keywords to try
	 *@return the function returns the maximum number of keywords that can be attempted; 
	 *this will be the same throughout the evaluation of the code, 
	 *so it may be more appropriate to have this as a "constant" attribute of some sort
	 *attemptNumber=-1 will call a special set of keywords for the monoatomic case     
	 */
	public int createGaussianPM3Input(String name, String directory, molFile p_molfile, int attemptNumber, String InChIaug, int multiplicity){
		//write a file with the input keywords
		QMInputWritable writer = new GaussianPM3InputWriter(name, directory, p_molfile, attemptNumber, multiplicity, InChIaug);

		File inputFile = writer.write();

		return QMInputWriter.maxAttemptNumber;
	}

	//creates MM4 input file and MM4 batch file in directory with filenames name.mm4 and name.com, respectively using MoleCoor
	//attemptNumber determines which keywords to try
	//the function returns the maximum number of keywords that can be attempted; this will be the same throughout the evaluation of the code, so it may be more appropriate to have this as a "constant" attribute of some sort
	public int createMM4Input(String name, String directory, molFile p_molfile, int attemptNumber, String InChIaug, int multiplicity){
		//write a file with the input keywords
		QMInputWritable writer = new MM4InputWriter(name, directory, p_molfile, attemptNumber, multiplicity, InChIaug);

		File inputFile = writer.write();

		return QMInputWriter.maxAttemptNumber;
	}

	
	/**
	 * creates MM4 rotor input file and MM4 batch file in directory 
	 * with filenames name.mm4roti and name.comi, respectively
	 * the function returns the set of rotor dihedral angles for 
	 * the minimum energy conformation 
	 */
	public double[] createMM4RotorInput(String name, String directory, ChemGraph p_chemgraph, int rotors){
		/*
		 * TODO should we include getter method of dihedral minima in interface QMInputWritable?
		 * don't think so...
		 * 
		 *  so for now, the implementation is used instead of the interface
		 */
		MM4HRInputWriter writer = new MM4HRInputWriter(name, directory, p_chemgraph, rotors);
		writer.write();
		return writer.getDihedralMinima();
	}

	//given x,y,z (cartesian) coordinates for dihedral1 atom, (rotor) atom 1, (rotor) atom 2, and dihedral2 atom, calculates the dihedral angle (in degrees, between +180 and -180) using the atan2 formula at http://en.wikipedia.org/w/index.php?title=Dihedral_angle&oldid=373614697
	public double calculateDihedral(double[] dihedral1, double[] atom1, double[] atom2, double[] dihedral2){
		//calculate the vectors between the atoms
		double[] b1 = {atom1[0]-dihedral1[0], atom1[1]-dihedral1[1], atom1[2]-dihedral1[2]};
		double[] b2 = {atom2[0]-atom1[0], atom2[1]-atom1[1], atom2[2]-atom1[2]};
		double[] b3 = {dihedral2[0]-atom2[0], dihedral2[1]-atom2[1], dihedral2[2]-atom2[2]};
		//calculate norm of b2
		double normb2 = Math.sqrt(b2[0]*b2[0]+b2[1]*b2[1]+b2[2]*b2[2]);
		//calculate necessary cross products
		double[] b1xb2 = {b1[1]*b2[2]-b1[2]*b2[1], b2[0]*b1[2]-b1[0]*b2[2], b1[0]*b2[1]-b1[1]*b2[0]};
		double[] b2xb3 = {b2[1]*b3[2]-b2[2]*b3[1], b3[0]*b2[2]-b2[0]*b3[2], b2[0]*b3[1]-b2[1]*b3[0]};
		//compute arguments for atan2 function (includes dot products)
		double y = normb2*(b1[0]*b2xb3[0]+b1[1]*b2xb3[1]+b1[2]*b2xb3[2]);//|b2|*b1.(b2xb3)
		double x = b1xb2[0]*b2xb3[0]+b1xb2[1]*b2xb3[1]+b1xb2[2]*b2xb3[2];//(b1xb2).(b2xb3)
		double dihedral = Math.atan2(y, x);
		//return dihedral*180.0/Math.PI;//converts from radians to degrees
		return Math.toDegrees(dihedral);//converts from radians to degrees
	}


	/**
	 * creates MOPAC PM3 input file in directory with filename name.mop by using OpenBabel to convert p_molfile
	 * attemptNumber determines which keywords to try
	 * the function returns the maximum number of keywords that can be attempted; this will be the same throughout the evaluation of the code, so it may be more appropriate to have this as a "constant" attribute of some sort
	 * unlike createGaussianPM3 input, this requires an additional input specifying the spin multiplicity (radical number + 1) for the species
	 * @param name
	 * @param directory
	 * @param p_molfile
	 * @param attemptNumber
	 * @param InChIaug
	 * @param multiplicity
	 * @return
	 */

	public int createMopacPM3Input(String name, String directory, molFile p_molfile, int attemptNumber, String InChIaug, int multiplicity){
		//write a file with the input keywords

		QMInputWritable writer = new MOPACPM3InputWriter(name, directory, p_molfile, attemptNumber, multiplicity);
		File inputFile = writer.write();

		return QMInputWriter.maxAttemptNumber;

	}


	/**
	 * name and directory are the name and directory for the input (and output) file;
	 * input is assumed to be preexisting and have the .gjf suffix
	 * returns an integer indicating success or failure of the Gaussian calculation: 1 for success, 0 for failure;
	 */
	public int runGaussian(String name, String directory){

		QMJobRunnable gaussianJob = new GaussianJob(name, directory);
		return gaussianJob.run();

	}

	//name and directory are the name and directory for the input (and output) file;
	//input script is assumed to be preexisting and have the .com suffix
	//returns an integer indicating success or failure of the calculation: 1 for success, 0 for failure
	public int runMM4(String name, String directory){

		QMJobRunnable mm4Job = new MM4Job(name, directory);
		return mm4Job.run();

	}

	//name and directory are the name and directory for the input (and output) file;
	//input script is assumed to be preexisting and have the .comi suffix where i is a number between 1 and rotors
	public void runMM4Rotor(String name, String directory, int rotors){
		for(int j=1;j<=rotors;j++){

			/**
			 * TODO for now, name for MM4HR job also contains some 
			 * input file extension and index j...
			 * 
			 * should be modified in the future!
			 */
			name = name+".com"+j;
			QMJobRunnable job = new MM4HRJob(name, directory);
			job.run();
		}
		return;
	}

	/**
	 * name and directory are the name and directory for the input (and output) file;
	 * input is assumed to be preexisting and have the .mop suffix
	 * returns an integer indicating success or failure of the MOPAC calculation: 1 for success, 0 for failure;
	 * this function is based on the Gaussian analogue
	 */
	public int runMOPAC(String name, String directory){

		QMJobRunnable job = new MOPACJob(name, directory);
		return job.run();

	}

	//parse the results using cclib and return a ThermoData object; name and directory indicate the location of the Gaussian .log file
	//may want to split this into several functions
	public ThermoData parseGaussianPM3(String name, String directory, ChemGraph p_chemGraph){

		QMParsable parser = new GaussianPM3Parser(name, directory, p_chemGraph);
		ThermoData result = parser.parse();
		result.setSource("Gaussian PM3 calculation");
		Logger.info("Thermo for " + name + ": "+ result.toString());//print result, at least for debugging purposes
		return result;
	}

	//parse the results using cclib and return a ThermoData object; name and directory indicate the location of the MM4 .mm4out file
	public ThermoData parseMM4(String name, String directory, ChemGraph p_chemGraph){

		QMParsable parser = new MM4Parser(name, directory, p_chemGraph);

		ThermoData result = parser.parse();

		result.setSource("MM4 calculation");
		Logger.info("Thermo for " + name + ": "+ result.toString());//print result, at least for debugging purposes
		return result;
	}

	//parse the results using cclib and CanTherm and return a ThermoData object; name and directory indicate the location of the MM4 .mm4out file
	//formerly known as parseMM4withForceMat
	public QMData performCanThermCalcs(String name, String directory, ChemGraph p_chemGraph, double[] dihedralMinima, boolean forceRRHO){
		//1. parse the MM4 file with cclib to get atomic number vector and geometry
		QMParser parser = new MM4Parser(name, directory, p_chemGraph, true);

		ThermoData data = parser.parse();
		
		IQMData qmdata = parser.getQMData();
		//unpack the needed results
		double energy = qmdata.energy;
		double stericEnergy = qmdata.stericEnergy;
		ArrayList freqs = qmdata.freqs;
		//2. compute H0/E0;  note that we will compute H0 for CanTherm by H0=Hf298(harmonicMM4)-(H298-H0)harmonicMM4, where harmonicMM4 values come from cclib parsing;  298.16 K is the standard temperature used by MM4; also note that Hthermal should include the R*T contribution (R*298.16 (enthalpy vs. energy difference)) and H0=E0 (T=0)	
		double T_MM4 = 298.16;
		energy *= QMConstants.Hartree_to_kcal;//convert from Hartree to kcal/mol
		stericEnergy *= QMConstants.Hartree_to_kcal;//convert from Hartree to kcal/mol
		double Hthermal = 5./2.*QMConstants.R*T_MM4/1000.;//enthalpy vs. energy difference(RT)+translation(3/2*R*T) contribution to thermal energy
		//rotational contribution
		if(p_chemGraph.getAtomNumber()==1) Hthermal += 0.0;
		else if (p_chemGraph.isLinear()) Hthermal += QMConstants.R*T_MM4/1000.;
		else Hthermal += 3./2.*QMConstants.R*T_MM4/1000;
		//vibrational contribution
		if(p_chemGraph.getAtomNumber()!=1)Hthermal+=QMConstants.R*calcVibH(freqs, T_MM4, QMConstants.h, QMConstants.k, QMConstants.c)/1000.;
		energy = energy - Hthermal;
		//3. write CanTherm input file
		//determine point group using the SYMMETRY Program
		String geom = qmdata.getSYMMETRYinput();
		String pointGroup = determinePointGroupUsingSYMMETRYProgram(geom);
		double sigmaCorr = getSigmaCorr(pointGroup);
		String canInp = "Calculation: Thermo\n";
		canInp += "Trange: 300 100 13\n";//temperatures from 300 to 1500 in increments of 100
		canInp += "Scale: 1.0\n";//scale factor of 1
		canInp += "Mol 1\n";
		if(p_chemGraph.getAtomNumber()==1) canInp += "ATOM\n";
		else if (p_chemGraph.isLinear()) canInp+="LINEAR\n";
		else canInp+="NONLINEAR\n";
		canInp += "GEOM MM4File " + name+".mm4out\n";//geometry file; ***special MM4 treatment in CanTherm; another option would be to use mm4opt file, but CanTherm code would need to be modified accordingly
		canInp += "FORCEC MM4File "+name+".fmat\n";//force constant file; ***special MM4 treatment in CanTherm
		if (forceRRHO) name = name + "_RRHO"; //"_RRHO" will be appended to InChIKey for RRHO calcs (though we want to use unmodified name in getQMDataWithCClib above, and in GEOM and FORCEC sections
		canInp += "ENERGY "+ energy +" MM4\n";//***special MM4 treatment in CanTherm
		canInp+="EXTSYM "+Math.exp(-sigmaCorr)+"\n";//***modified treatment in CanTherm; traditional EXTSYM integer replaced by EXTSYM double, to allow fractional values that take chirality into account
		canInp+="NELEC 1\n";//multiplicity = 1; all cases we consider will be non-radicals
		int rotors = p_chemGraph.getInternalRotor();
		String rotInput = null;
		if(!useHindRot || rotors==0 || forceRRHO) canInp += "ROTORS 0\n";//do not consider hindered rotors
		else{
			int rotorCount = 0;
			canInp+="ROTORS "+rotors+ " "+name+".rotinfo\n";
			canInp+="POTENTIAL separable mm4files_inertia\n";//***special MM4 treatment in canTherm;
			String rotNumbersLine=""+stericEnergy;//the line that will contain V0 (kcal/mol), and all the dihedral minima (degrees)
			rotInput = "L1: 1 2 3\n";
			LinkedHashMap rotorInfo = p_chemGraph.getInternalRotorInformation();
			Iterator iter = rotorInfo.keySet().iterator();
			while(iter.hasNext()){
				rotorCount++;
				int[] rotorAtoms = (int[])iter.next();
				LinkedHashSet rotatingGroup = (LinkedHashSet)rotorInfo.get(rotorAtoms);
				Iterator iterGroup = rotatingGroup.iterator();
				rotInput += "L2: " + rotorAtoms[4] +" "+rotorAtoms[1]+ " "+ rotorAtoms[2];//note: rotorAtoms[4] is the rotorSymmetry number as estimated by calculateRotorSymmetryNumber; this will be taken into account elsewhere
				while (iterGroup.hasNext()){//print the atoms associated with rotorAtom 2
					Integer id = (Integer)iterGroup.next();
					if(id != rotorAtoms[2]) rotInput+= " "+id;//don't print atom2
				}
				rotInput += "\n";
				canInp+=name + ".mm4rotopt"+rotorCount+" ";// potential files will be named as name.mm4rotopti
				rotNumbersLine+=" "+dihedralMinima[rotorCount-1];
			}
			canInp+="\n"+rotNumbersLine+"\n";
		}
		canInp+="0\n";//no bond-additivity corrections
		try{
			File canFile=new File(directory+"/"+name+".can");
			FileWriter fw = new FileWriter(canFile);
			fw.write(canInp);
			fw.close();
			if(rotInput !=null){//write the rotor information
				File rotFile=new File(directory+"/"+name+".rotinfo");
				FileWriter fwr = new FileWriter(rotFile);
				fwr.write(rotInput);
				fwr.close();
			}
		}
		catch(Exception e){
			String err = "Error in writing CanTherm input \n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}
		//4 call CanTherm 
		QMJobRunnable canthermJob = new CanThermJob(name, directory);
		canthermJob.run();
		
		return qmdata;
	}

	public ThermoData parseCanThermFile(String name, String directory, ChemGraph p_chemGraph){
		//5. read CanTherm output
		Double Hf298 = null;
		Double S298 = null;
		Double Cp300 = null;
		Double Cp400 = null;
		Double Cp500 = null;
		Double Cp600 = null;
		Double Cp800 = null;
		Double Cp1000 = null;
		Double Cp1500 = null;
		File file = new File(directory+"/"+name+".canout");
		try{
			FileReader in = new FileReader(file);
			BufferedReader reader = new BufferedReader(in);
			String line=reader.readLine();
			while(!line.startsWith("Hf298 S298 Cps:")){//get to the end of the file with the data we want
				line=reader.readLine();
			}
			String[] split = reader.readLine().trim().split("\\s+");//read the next line, which contains the information we want
			Hf298 = Double.parseDouble(split[0]);
			S298 = Double.parseDouble(split[1]);
			Cp300 = Double.parseDouble(split[2]);
			Cp400 = Double.parseDouble(split[3]);
			Cp500 = Double.parseDouble(split[4]);
			Cp600 = Double.parseDouble(split[5]);
			Cp800 = Double.parseDouble(split[7]);
			Cp1000 = Double.parseDouble(split[9]);
			Cp1500 = Double.parseDouble(split[14]);
			reader.close();
			in.close();
		}
		catch(Exception e){
			String err = "Error in reading CanTherm .canout file \n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}

		ThermoData result = new ThermoData(Hf298,S298,Cp300,Cp400,Cp500,Cp600,Cp800,Cp1000,Cp1500,3,1,1,"MM4 calculation; includes CanTherm analysis of force-constant matrix");//this includes rough estimates of uncertainty
		result.setSource("MM4 calculation with CanTherm analysis");
		Logger.info("Thermo for " + name + ": "+ result.toString());//print result, at least for debugging purposes

		return result;
	}

	//parse the results using cclib and return a ThermoData object; name and directory indicate the location of the MOPAC .out file
	public ThermoData parseMopacPM3(String name, String directory, ChemGraph p_chemGraph){

		QMParsable parser = new MOPACPM3Parser(name, directory, p_chemGraph);

		ThermoData result = parser.parse();
		result.setSource("MOPAC PM3 calculation");
		Logger.info("Thermo for " + name + ": "+ result.toString());//print result, at least for debugging purposes
		return result;
	}

	//gets the statistical correction for S in dimensionless units (divided by R)
	public double getSigmaCorr(String pointGroup){
		double sigmaCorr=0;
		//determine statistical correction factor for 1. external rotational symmetry (affects rotational partition function) and 2. chirality (will add R*ln2 to entropy) based on point group
		//ref: http://cccbdb.nist.gov/thermo.asp
		//assumptions below for Sn, T, Th, O, I seem to be in line with expectations based on order reported at: http://en.wikipedia.org/w/index.php?title=List_of_character_tables_for_chemically_important_3D_point_groups&oldid=287261611 (assuming order = symmetry number * 2 (/2 if chiral))...this appears to be true for all point groups I "know" to be correct
		//minor concern: does SYMMETRY appropriately calculate all Sn groups considering 2007 discovery of previous errors in character tables (cf. Wikipedia article above)
		if (pointGroup.equals("C1")) sigmaCorr=+Math.log(2.);//rot. sym. = 1, chiral
		else if (pointGroup.equals("Cs")) sigmaCorr=0; //rot. sym. = 1
		else if (pointGroup.equals("Ci")) sigmaCorr=0; //rot. sym. = 1
		else if (pointGroup.equals("C2")) sigmaCorr=0;//rot. sym. = 2, chiral (corrections cancel)
		else if (pointGroup.equals("C3")) sigmaCorr=+Math.log(2.)-Math.log(3.);//rot. sym. = 3, chiral
		else if (pointGroup.equals("C4")) sigmaCorr=+Math.log(2.)-Math.log(4.);//rot. sym. = 4, chiral
		else if (pointGroup.equals("C5")) sigmaCorr=+Math.log(2.)-Math.log(5.);//rot. sym. = 5, chiral
		else if (pointGroup.equals("C6")) sigmaCorr=+Math.log(2.)-Math.log(6.);//rot. sym. = 6, chiral
		else if (pointGroup.equals("C7")) sigmaCorr=+Math.log(2.)-Math.log(7.);//rot. sym. = 7, chiral
		else if (pointGroup.equals("C8")) sigmaCorr=+Math.log(2.)-Math.log(8.);//rot. sym. = 8, chiral
		else if (pointGroup.equals("D2")) sigmaCorr=+Math.log(2.)-Math.log(4.);//rot. sym. = 4, chiral
		else if (pointGroup.equals("D3")) sigmaCorr=+Math.log(2.)-Math.log(6.);//rot. sym. = 6, chiral
		else if (pointGroup.equals("D4")) sigmaCorr=+Math.log(2.)-Math.log(8.);//rot. sym. = 8, chiral
		else if (pointGroup.equals("D5")) sigmaCorr=+Math.log(2.)-Math.log(10.);//rot. sym. = 10, chiral
		else if (pointGroup.equals("D6")) sigmaCorr=+Math.log(2.)-Math.log(12.);//rot. sym. = 12, chiral
		else if (pointGroup.equals("D7")) sigmaCorr=+Math.log(2.)-Math.log(14.);//rot. sym. = 14, chiral
		else if (pointGroup.equals("D8")) sigmaCorr=+Math.log(2.)-Math.log(16.);//rot. sym. = 16, chiral
		else if (pointGroup.equals("C2v")) sigmaCorr=-Math.log(2.);//rot. sym. = 2
		else if (pointGroup.equals("C3v")) sigmaCorr=-Math.log(3.);//rot. sym. = 3
		else if (pointGroup.equals("C4v")) sigmaCorr=-Math.log(4.);//rot. sym. = 4
		else if (pointGroup.equals("C5v")) sigmaCorr=-Math.log(5.);//rot. sym. = 5
		else if (pointGroup.equals("C6v")) sigmaCorr=-Math.log(6.);//rot. sym. = 6
		else if (pointGroup.equals("C7v")) sigmaCorr=-Math.log(7.);//rot. sym. = 7
		else if (pointGroup.equals("C8v")) sigmaCorr=-Math.log(8.);//rot. sym. = 8
		else if (pointGroup.equals("C2h")) sigmaCorr=-Math.log(2.);//rot. sym. = 2
		else if (pointGroup.equals("C3h")) sigmaCorr=-Math.log(3.);//rot. sym. = 3
		else if (pointGroup.equals("C4h")) sigmaCorr=-Math.log(4.);//rot. sym. = 4
		else if (pointGroup.equals("C5h")) sigmaCorr=-Math.log(5.);//rot. sym. = 5
		else if (pointGroup.equals("C6h")) sigmaCorr=-Math.log(6.);//rot. sym. = 6
		else if (pointGroup.equals("C7h")) sigmaCorr=-Math.log(7.);//rot. sym. = 7
		else if (pointGroup.equals("C8h")) sigmaCorr=-Math.log(8.);//rot. sym. = 8
		else if (pointGroup.equals("D2h")) sigmaCorr=-Math.log(4.);//rot. sym. = 4
		else if (pointGroup.equals("D3h")) sigmaCorr=-Math.log(6.);//rot. sym. = 6
		else if (pointGroup.equals("D4h")) sigmaCorr=-Math.log(8.);//rot. sym. = 8
		else if (pointGroup.equals("D5h")) sigmaCorr=-Math.log(10.);//rot. sym. = 10
		else if (pointGroup.equals("D6h")) sigmaCorr=-Math.log(12.);//rot. sym. = 12
		else if (pointGroup.equals("D7h")) sigmaCorr=-Math.log(14.);//rot. sym. = 14
		else if (pointGroup.equals("D8h")) sigmaCorr=-Math.log(16.);//rot. sym. = 16
		else if (pointGroup.equals("D2d")) sigmaCorr=-Math.log(4.);//rot. sym. = 4
		else if (pointGroup.equals("D3d")) sigmaCorr=-Math.log(6.);//rot. sym. = 6
		else if (pointGroup.equals("D4d")) sigmaCorr=-Math.log(8.);//rot. sym. = 8
		else if (pointGroup.equals("D5d")) sigmaCorr=-Math.log(10.);//rot. sym. = 10
		else if (pointGroup.equals("D6d")) sigmaCorr=-Math.log(12.);//rot. sym. = 12
		else if (pointGroup.equals("D7d")) sigmaCorr=-Math.log(14.);//rot. sym. = 14
		else if (pointGroup.equals("D8d")) sigmaCorr=-Math.log(16.);//rot. sym. = 16
		else if (pointGroup.equals("S4")) sigmaCorr=-Math.log(2.);//rot. sym. = 2 ;*** assumed achiral
		else if (pointGroup.equals("S6")) sigmaCorr=-Math.log(3.);//rot. sym. = 3 ;*** assumed achiral
		else if (pointGroup.equals("S8")) sigmaCorr=-Math.log(4.);//rot. sym. = 4 ;*** assumed achiral
		else if (pointGroup.equals("T")) sigmaCorr=+Math.log(2.)-Math.log(12.);//rot. sym. = 12, *** assumed chiral
		else if (pointGroup.equals("Th")) sigmaCorr=-Math.log(12.);//***assumed rot. sym. = 12
		else if (pointGroup.equals("Td")) sigmaCorr=-Math.log(12.);//rot. sym. = 12
		else if (pointGroup.equals("O")) sigmaCorr=+Math.log(2.)-Math.log(24.);//***assumed rot. sym. = 24, chiral
		else if (pointGroup.equals("Oh")) sigmaCorr=-Math.log(24.);//rot. sym. = 24
		else if (pointGroup.equals("Cinfv")) sigmaCorr=0;//rot. sym. = 1
		else if (pointGroup.equals("Dinfh")) sigmaCorr=-Math.log(2.);//rot. sym. = 2
		else if (pointGroup.equals("I")) sigmaCorr=+Math.log(2.)-Math.log(60.);//***assumed rot. sym. = 60, chiral
		else if (pointGroup.equals("Ih")) sigmaCorr=-Math.log(60.);//rot. sym. = 60
		else if (pointGroup.equals("Kh")) sigmaCorr=0;//arbitrarily set to zero...one could argue that it is infinite; apparently this is the point group of a single atom (cf. http://www.cobalt.chem.ucalgary.ca/ps/symmetry/tests/G_Kh); this should not have a rotational partition function, and we should not use the symmetry correction in this case
		else{//this point should not be reached, based on checks performed in determinePointGroupUsingSYMMETRYProgram
			Logger.critical("Unrecognized point group: "+ pointGroup);
			System.exit(0);
		}

		return sigmaCorr;
	}


	//determine the point group using the SYMMETRY program (http://www.cobalt.chem.ucalgary.ca/ps/symmetry/)
	//required input is a line with number of atoms followed by lines for each atom including atom number and x,y,z coordinates
	//finalTol determines how loose the point group criteria are; values are comparable to those specifed in the GaussView point group interface
	//public String determinePointGroupUsingSYMMETRYProgram(String geom, double finalTol){
	public String determinePointGroupUsingSYMMETRYProgram(String geom){
		int attemptNumber = 1;
		int maxAttemptNumber = 4;
		boolean pointGroupFound=false;
		//write the input file
		try {
			File inputFile=new File(qmfolder+"symminput.txt");//SYMMETRY program directory
			FileWriter fw = new FileWriter(inputFile);
			fw.write(geom);
			fw.close();
		} catch (IOException e) {
			String err = "Error writing input file for point group calculation";
			err += e.toString();
			Logger.critical(err);
			System.exit(0);
		}
		String result = "";
		String command = "";
		while (attemptNumber<=maxAttemptNumber && !pointGroupFound){
			//call the program and read the result
			result = "";
			String [] lineArray;
			try{ 
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
				Process symmProc = Runtime.getRuntime().exec(command);
				//check for errors and display the error if there is one
				InputStream is = symmProc.getInputStream();
				InputStreamReader isr = new InputStreamReader(is);
				BufferedReader br = new BufferedReader(isr);
				String line=null;
				while ( (line = br.readLine()) != null) {
					if(line.startsWith("It seems to be the ")){//last line, ("It seems to be the [x] point group") indicates point group
						lineArray = line.split(" ");//split the line around spaces
						result = lineArray[5];//point group string should be the 6th word
					}
				}
				int exitValue = symmProc.waitFor();
				symmProc.getErrorStream().close();
				symmProc.getOutputStream().close();
				br.close();
				isr.close();
				is.close();
			}
			catch(Exception e){
				String err = "Error in running point group calculation process using SYMMETRY \n";
				err += e.toString();
				Logger.logStackTrace(e);
				System.exit(0);
			}
			//check for a recognized point group
			if (result.equals("C1")||result.equals("Cs")||result.equals("Ci")||result.equals("C2")||result.equals("C3")||result.equals("C4")||result.equals("C5")||result.equals("C6")||result.equals("C7")||result.equals("C8")||result.equals("D2")||result.equals("D3")||result.equals("D4")||result.equals("D5")||result.equals("D6")||result.equals("D7")||result.equals("D8")||result.equals("C2v")||result.equals("C3v")||result.equals("C4v")||result.equals("C5v")||result.equals("C6v")||result.equals("C7v")||result.equals("C8v")||result.equals("C2h")||result.equals("C3h")||result.equals("C4h")||result.equals("C5h")||result.equals("C6h")||result.equals("C7h")||result.equals("C8h")||result.equals("D2h")||result.equals("D3h")||result.equals("D4h")||result.equals("D5h")||result.equals("D6h")||result.equals("D7h")||result.equals("D8h")||result.equals("D2d")||result.equals("D3d")||result.equals("D4d")||result.equals("D5d")||result.equals("D6d")||result.equals("D7d")||result.equals("D8d")||result.equals("S4")||result.equals("S6")||result.equals("S8")||result.equals("T")||result.equals("Th")||result.equals("Td")||result.equals("O")||result.equals("Oh")||result.equals("Cinfv")||result.equals("Dinfh")||result.equals("I")||result.equals("Ih")||result.equals("Kh")) pointGroupFound=true;
			else{
				if(attemptNumber < maxAttemptNumber) Logger.info("Attempt number "+attemptNumber+" did not identify a recognized point group (" +result+"). Will retry with looser point group criteria.");
				else{
					Logger.critical("Final attempt number "+attemptNumber+" did not identify a recognized point group (" +result+"). Exiting.");
					System.exit(0);
				}
				attemptNumber++;
			}
		} 
		Logger.info("Point group: "+ result);//print result, at least for debugging purposes

		return result;
	}

	//gmagoon 6/23/10
	//calculate the vibrational contribution (divided by R, units of K) at temperature, T, in Kelvin to Hthermal (ZPE not included)
	//p_freqs in cm^-1; c in cm/s; k in J/K; h in J-s
	//we need to ignore zero frequencies, as MM4 does, to avoid dividing by zero; based on equation for Ev in http://www.gaussian.com/g_whitepap/thermo.htm; however, we ignore zero point contribution to be consistent with the input that CanTherm currently takes; note, however, that this could introduce some small inaccuracies, as the frequencies may be slightly different in CanTherm vs. MM4, particularly for a frequency < 7.7 cm^-1 (treated as zero in MM4)
	public double calcVibH(ArrayList p_freqs, double p_T, double h, double k, double c){
		double Hcontrib = 0;
		double dr;
		for(int i=0; i < p_freqs.size(); i++){
			double freq = (Double)p_freqs.get(i);
			if(freq > 0.0){//ignore zero frequencies as MM4 does
				dr = h*c*freq/(k*p_T); //frequently used dimensionless ratio
				Hcontrib = Hcontrib + dr*p_T/(Math.exp(dr) - 1.);
			}
		}

		return Hcontrib;
	}


	//determine the QM filename (element 0) and augmented InChI (element 1) for a ChemGraph
	//QM filename is InChIKey appended with mult3, mult4, mult5, or mult6 for multiplicities of 3 or higher
	//augmented InChI is InChI appended with /mult3, /mult4, /mult5, or /mult6 for multiplicities of 3 or higher
	public String [] getQMFileName(ChemGraph p_chemGraph){
		String [] result = new String[2];
		result[0] = p_chemGraph.getModifiedInChIKeyAnew();//need to generate InChI and key anew because ChemGraph may have changed (in particular, adding/removing hydrogens in HBI process)
		result[1] = p_chemGraph.getModifiedInChIAnew();

		return result;
	}

	//returns true if a Gaussian file for the given name and directory (.log suffix) exists and indicates successful completion (same criteria as used after calculation runs); terminates if the InChI doesn't match the InChI in the file or if there is no InChI in the file; returns false otherwise
	public boolean successfulGaussianResultExistsQ(String name, String directory, String InChIaug){
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
	public boolean successfulMopacResultExistsQ(String name, String directory, String InChIaug){
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
	public boolean successfulMM4ResultExistsQ(String name, String directory, String InChIaug){
		//part of the code is taken from analogous code for MOPAC (first ~half) and Gaussian (second ~half)
		//look in the output file to check for the successful termination of the calculation (assumed to be successful if "description of vibrations appears)
		int failureFlag=1;//flag (1 or 0) indicating whether the MM4 job failed
		int failureOverrideFlag=0;//flag (1 or 0) to override success as measured by failureFlag
		File file = new File(directory+"/"+name+".mm4out");
		File canFile = new File(directory+"/"+name+".canout");
		int InChIMatch=0;//flag (1 or 0) indicating whether the InChI in the file matches InChIaug; this can only be 1 if InChIFound is also 1;
		int InChIFound=0;//flag (1 or 0) indicating whether an InChI was found in the log file
		int InChIPartialMatch=0;//flag (1 or 0) indicating whether the InChI in the log file is a substring of the InChI in RMG's memory
		if(useCanTherm){//if we are using CanTherm, check whether a CanTherm output file exists...if it does, we will continue on, otherwise, we will rerun calculations (including MM4 calculation) from scratch to ensure our atom numbering is consistent; note: if the .canout file exists, we still want to check for the actual MM4 file even if we are using CanTherm and reading CanTherm output because 1. it ensures we have the correct species and don't have an InChI collision 2. it is needed for getting the geometry which is needed for the symmetry number corrections applied to CanTherm output (which doesn't include symmetry number considerations)
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
						if (useCanTherm){//zero frequencies are only acceptable when CanTherm is used
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

	//    //checks the MOPAC file for negative frequencies
	//    public boolean MopacFileContainsNegativeFreqsQ(String name, String directory){
	//        boolean negativeFreq=false;
	//        
	//        //code below copied from parseMopacPM3()
	//        String command = "c:/Python25/python.exe c:/Python25/MopacPM3ParsingScript.py ";//this should eventually be modified for added generality
	//        String logfilepath=directory+"/"+name+".out";
	//        command=command.concat(logfilepath);
	//        
	//        //much of code below is copied from calculateThermoFromPM3Calc()
	//        //parse the Mopac file using cclib
	//        int natoms = 0; //number of atoms from Mopac file; in principle, this should agree with number of chemGraph atoms
	//        ArrayList atomicNumber = new ArrayList(); //vector of atomic numbers (integers) (apparently Vector is thread-safe; cf. http://answers.yahoo.com/question/index?qid=20081214065127AArZDT3; ...should I be using this instead?)
	//        ArrayList x_coor = new ArrayList(); //vectors of x-, y-, and z-coordinates (doubles) (Angstroms) (in order corresponding to above atomic numbers)
	//        ArrayList y_coor = new ArrayList();
	//        ArrayList z_coor = new ArrayList();
	//        double energy = 0; //PM3 energy (Hf298) in Hartree (***note: in the case of MOPAC, the MOPAC file will contain in units of kcal/mol, but modified ccLib will return in Hartree)
	//        double molmass = 0; //molecular mass in amu
	//        ArrayList freqs = new ArrayList(); //list of frequencies in units of cm^-1
	//        double rotCons_1 = 0;//rotational constants in (1/s)
	//        double rotCons_2 = 0;
	//        double rotCons_3 = 0; 
	//        //int gdStateDegen = p_chemGraph.getRadicalNumber()+1;//calculate ground state degeneracy from the number of radicals; this should give the same result as spin multiplicity in Gaussian input file (and output file), but we do not explicitly check this (we could use "mult" which cclib reads in if we wanted to do so); also, note that this is not always correct, as there can apparently be additional spatial degeneracy for non-symmetric linear molecules like OH radical (cf. http://cccbdb.nist.gov/thermo.asp)
	//        try{   
	//            File runningdir=new File(directory);
	//            Process cclibProc = Runtime.getRuntime().exec(command, null, runningdir);
	//            //read the stdout of the process, which should contain the desired information in a particular format
	//            InputStream is = cclibProc.getInputStream();
	//            InputStreamReader isr = new InputStreamReader(is);
	//            BufferedReader br = new BufferedReader(isr);
	//            String line=null;
	//            //example output:
	////            C:\Python25>python.exe GaussianPM3ParsingScript.py TEOS.out
	////            33
	////            [ 6  6  8 14  8  6  6  8  6  6  8  6  6  1  1  1  1  1  1  1  1  1  1  1  1
	////              1  1  1  1  1  1  1  1]
	////            [[ 2.049061 -0.210375  3.133106]
	////             [ 1.654646  0.321749  1.762752]
	////             [ 0.359284 -0.110429  1.471465]
	////             [-0.201871 -0.013365 -0.12819 ]
	////             [ 0.086307  1.504918 -0.82893 ]
	////             [-0.559186  2.619928 -0.284003]
	////             [-0.180246  3.839463 -1.113029]
	////             [ 0.523347 -1.188305 -1.112765]
	////             [ 1.857584 -1.018167 -1.495088]
	////             [ 2.375559 -2.344392 -2.033403]
	////             [-1.870397 -0.297297 -0.075427]
	////             [-2.313824 -1.571765  0.300245]
	////             [-3.83427  -1.535927  0.372171]
	////             [ 1.360346  0.128852  3.917699]
	////             [ 2.053945 -1.307678  3.160474]
	////             [ 3.055397  0.133647  3.403037]
	////             [ 1.677262  1.430072  1.750899]
	////             [ 2.372265 -0.029237  0.985204]
	////             [-0.245956  2.754188  0.771433]
	////             [-1.656897  2.472855 -0.287156]
	////             [-0.664186  4.739148 -0.712606]
	////             [-0.489413  3.734366 -2.161038]
	////             [ 0.903055  4.016867 -1.112198]
	////             [ 1.919521 -0.229395 -2.269681]
	////             [ 2.474031 -0.680069 -0.629949]
	////             [ 2.344478 -3.136247 -1.273862]
	////             [ 1.786854 -2.695974 -2.890647]
	////             [ 3.41648  -2.242409 -2.365094]
	////             [-1.884889 -1.858617  1.28054 ]
	////             [-1.976206 -2.322432 -0.440995]
	////             [-4.284706 -1.26469  -0.591463]
	////             [-4.225999 -2.520759  0.656131]
	////             [-4.193468 -0.809557  1.112677]]
	////            -14.1664924726
	////            [    9.9615    18.102     27.0569    31.8459    39.0096    55.0091
	////                66.4992    80.4552    86.4912   123.3551   141.6058   155.5448
	////               159.4747   167.0013   178.5676   207.3738   237.3201   255.3487
	////               264.5649   292.867    309.4248   344.6503   434.8231   470.2074
	////               488.9717   749.1722   834.257    834.6594   837.7292   839.6352
	////               887.9767   892.9538   899.5374   992.1851  1020.6164  1020.8671
	////              1028.3897  1046.7945  1049.1768  1059.4704  1065.1505  1107.4001
	////              1108.1567  1109.0466  1112.6677  1122.7785  1124.4315  1128.4163
	////              1153.3438  1167.6705  1170.9627  1174.9613  1232.1826  1331.8459
	////              1335.3932  1335.8677  1343.9556  1371.37    1372.8127  1375.5428
	////              1396.0344  1402.4082  1402.7554  1403.2463  1403.396   1411.6946
	////              1412.2456  1412.3519  1414.5982  1415.3613  1415.5698  1415.7993
	////              1418.5409  2870.7446  2905.3132  2907.0361  2914.1662  2949.2646
	////              2965.825   2967.7667  2971.5223  3086.3849  3086.3878  3086.6448
	////              3086.687   3089.2274  3089.4105  3089.4743  3089.5841  3186.0753
	////              3186.1375  3186.3511  3186.365 ]
	////            [ 0.52729  0.49992  0.42466]
	////note: above example has since been updated to print molecular mass; also frequency and atomic number format has been updated
	//            String [] stringArray;
	//            natoms = Integer.parseInt(br.readLine());//read line 1: number of atoms
	//            stringArray = br.readLine().replace("[", "").replace("]","").trim().split(",\\s+");//read line 2: the atomic numbers (first removing braces)
	//           // line = br.readLine().replace("[", "").replace("]","");//read line 2: the atomic numbers (first removing braces)
	//           // StringTokenizer st = new StringTokenizer(line); //apprently the stringTokenizer class is deprecated, but I am having trouble getting the regular expressions to work properly
	//            for(int i=0; i < natoms; i++){
	//               // atomicNumber.add(i,Integer.parseInt(stringArray[i]));
	//                atomicNumber.add(i,Integer.parseInt(stringArray[i]));
	//            }
	//            for(int i=0; i < natoms; i++){
	//                stringArray = br.readLine().replace("[", "").replace("]","").trim().split("\\s+");//read line 3+: coordinates for atom i; used /s+ for split; using spaces with default limit of 0 was giving empty string
	//                x_coor.add(i,Double.parseDouble(stringArray[0]));
	//                y_coor.add(i,Double.parseDouble(stringArray[1]));
	//                z_coor.add(i,Double.parseDouble(stringArray[2]));
	//            }
	//            energy = Double.parseDouble(br.readLine());//read next line: energy
	//            molmass = Double.parseDouble(br.readLine());//read next line: molecular mass (in amu)
	//            if (natoms>1){//read additional info for non-monoatomic species
	//                stringArray = br.readLine().replace("[", "").replace("]","").trim().split(",\\s+");//read next line: frequencies
	//                for(int i=0; i < stringArray.length; i++){
	//                    freqs.add(i,Double.parseDouble(stringArray[i]));
	//                }
	//                stringArray = br.readLine().replace("[", "").replace("]","").trim().split("\\s+");//read next line rotational constants (converting from GHz to Hz in the process)
	//                rotCons_1 = Double.parseDouble(stringArray[0])*1000000000;
	//                rotCons_2 = Double.parseDouble(stringArray[1])*1000000000;
	//                rotCons_3 = Double.parseDouble(stringArray[2])*1000000000;
	//            }
	//            while ( (line = br.readLine()) != null) {
	//                //do nothing (there shouldn't be any more information, but this is included to get all the output)
	//            }
	//            int exitValue = cclibProc.waitFor();
	//        }
	//        catch (Exception e) {
	//				Logger.logStackTrace(e);
	//            String err = "Error in running ccLib Python process \n";
	//            err += e.toString();
	//            Logger.critical(err);
	//            System.exit(0);
	//        } 
	//       
	//        //start of code "new" to this function (aside from initialization of negativeFreq)
	//        if(natoms > 0){
	//            for (int i=0; i<freqs.size(); i++){
	//                if((Double)freqs.get(i) < 0) negativeFreq = true;
	//            }
	//        }
	//        return negativeFreq;
	//    }
	//## operation initGAGroupLibrary()
	protected void initGAGroupLibrary() {
		//#[ operation initGAGroupLibrary()
		thermoLibrary = ThermoGAGroupLibrary.getINSTANCE();
		//#]
	}


}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\QMTP.java
 *********************************************************************/
