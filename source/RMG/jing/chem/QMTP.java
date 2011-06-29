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

import jing.chem.qmtp.CanThermInputWriter;
import jing.chem.qmtp.CanThermJob;
import jing.chem.qmtp.CanThermParser;
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
import jing.chem.qmtp.ThreeDMolFileCreator;
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
				ThreeDMolFileCreator creator = new ThreeDMolFileCreator(name, p_chemGraph);
				molFile p_3dfile = creator.create();
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
						/*
						 * name and directory are the name and directory for the input (and output) file;
						 * input is assumed to be preexisting and have the .gjf suffix
						 * returns an integer indicating success or failure of the Gaussian calculation: 1 for success, 0 for failure;
						 */
						QMJobRunnable gaussianJob = new GaussianJob(name, directory);
						successFlag = gaussianJob.run();
					}
					else if (qmProgram.equals("mopac") || qmProgram.equals("both")){
						maxAttemptNumber = createMopacPM3Input(name, directory, p_3dfile, attemptNumber, InChIaug, multiplicity);
						/*
						 * name and directory are the name and directory for the input (and output) file;
						 * input is assumed to be preexisting and have the .mop suffix
						 * returns an integer indicating success or failure of the MOPAC calculation: 1 for success, 0 for failure;
						 * this function is based on the Gaussian analogue
						 */
						QMJobRunnable job = new MOPACJob(name, directory); 
						successFlag = job.run();
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
				//parse the results using cclib and return a ThermoData object; name and directory indicate the location of the Gaussian .log file
				//may want to split this into several functions
				QMParsable parser = new GaussianPM3Parser(name, directory, p_chemGraph);
				result = parser.parse();
				result.setSource("Gaussian PM3 calculation");
				Logger.info("Thermo for " + name + ": "+ result.toString());//print result, at least for debugging purposes
				
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
				ThreeDMolFileCreator creator = new ThreeDMolFileCreator(name, p_chemGraph);
				molFile p_3dfile = creator.create();
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
					//name and directory are the name and directory for the input (and output) file;
					//input script is assumed to be preexisting and have the .com suffix
					//returns an integer indicating success or failure of the calculation: 1 for success, 0 for failure
					QMJobRunnable mm4Job = new MM4Job(name, directory);
					successFlag = mm4Job.run();
					//new IF block to check success
					if(successFlag==1){
						Logger.info("Attempt #"+attemptNumber + " on species " + name + " ("+InChIaug+") succeeded.");
						//run rotor calculations if necessary
						int rotors = p_chemGraph.getInternalRotor();
						if(useHindRot && rotors > 0){
							//we should re-run scans even if pre-existing scans exist because atom numbering may change from case to case; a better solution would be to check for stored CanTherm output and use that if available
							Logger.info("Running rotor scans on "+name+"...");
							dihedralMinima = createMM4RotorInput(name, directory, p_chemGraph, rotors);//we don't worry about checking InChI here; if there is an InChI mismatch it should be caught
							//name and directory are the name and directory for the input (and output) file;
							//input script is assumed to be preexisting and have the .comi suffix where i is a number between 1 and rotors
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
				CanThermParser canthermParser = new CanThermParser(name, directory, p_chemGraph);
				result = canthermParser.parse();
				if (p_chemGraph.getInternalRotor()>0) {
					//print the RRHO result for comparison
					canthermParser = new CanThermParser(name+"_RRHO", directory, p_chemGraph);
				}
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
	public IQMData performCanThermCalcs(String name, String directory, ChemGraph p_chemGraph, double[] dihedralMinima, boolean forceRRHO){
		// parse the MM4 file with cclib to get atomic number vector and geometry
		QMParser parser = new MM4Parser(name, directory, p_chemGraph, true);

		/*
		 * TODO if we don't call the parse() method the IQMData object
		 * would not be created...
		 */
		ThermoData data = parser.parse();
		
		IQMData qmdata = parser.getQMData();
		//unpack the needed results and write cantherm input file
		QMInputWritable canthermWriter = new CanThermInputWriter(name, directory, p_chemGraph, qmdata, dihedralMinima, forceRRHO);
		File canThermFile = canthermWriter.write();
		
		// call CanTherm 
		QMJobRunnable canthermJob = new CanThermJob(name, directory);
		canthermJob.run();
		
		return qmdata;
	}

	//parse the results using cclib and return a ThermoData object; name and directory indicate the location of the MOPAC .out file
	public ThermoData parseMopacPM3(String name, String directory, ChemGraph p_chemGraph){

		QMParsable parser = new MOPACPM3Parser(name, directory, p_chemGraph);

		ThermoData result = parser.parse();
		result.setSource("MOPAC PM3 calculation");
		Logger.info("Thermo for " + name + ": "+ result.toString());//print result, at least for debugging purposes
		return result;
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
