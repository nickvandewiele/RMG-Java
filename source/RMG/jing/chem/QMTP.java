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
import jing.chem.qmtp.QMJobRunnable;
import jing.chem.qmtp.QMParsable;
import jing.chem.qmtp.QMParser;
import jing.chem.qmtp.QMVerifier;
import jing.chem.qmtp.ThreeDMolFileCreator;
import jing.chemUtil.*;
import jing.param.*;

import java.io.File;
import jing.rxnSys.Logger;

//quantum mechanics thermo property estimator; analog of GATP
public class QMTP implements GeneralGAPP {

	private static QMTP INSTANCE = new QMTP();		//## attribute INSTANCE
	protected static PrimaryThermoLibrary primaryLibrary;//Note: may be able to separate this out into GeneralGAPP, as this is common to both GATP and QMTP
	public static String qmfolder= "QMfiles/";
	//   protected static HashMap library;		//as above, may be able to move this and associated functions to GeneralGAPP (and possibly change from "x implements y" to "x extends y"), as it is common to both GATP and QMTP
	protected ThermoGAGroupLibrary thermoLibrary; //needed for HBI

	/*
	 * the qmprogram can be 
	 * "mopac", 
	 * "gaussian03",
	 *  "both" (MOPAC and Gaussian), 
	 *  "mm4",
	 *  "mm4hr"
	 */
	public static String qmprogram= "both";

	/*
	 * the computational method called by one of the QM Programs.
	 * can be:
	 * "PM3",
	 * "MM4",
	 * "MM4HR"
	 * 
	 */
	public static String qmMethod;

	public Map<String, Integer> mapMaxAttemptNumber;
	public static boolean usePolar = false; //use polar keyword in MOPAC
	public static boolean useCanTherm = true; //whether to use CanTherm in MM4 cases for interpreting output via force-constant matrix; this will hopefully avoid zero frequency issues
	public static boolean useHindRot = false;//whether to use HinderedRotor scans with MM4 (requires useCanTherm=true)
	//## operation QMTP()
	private QMTP() {
		/*
		 * set global attribute qmMethod
		 * 
		 * Right now, not many choices are available.
		 */
		if(qmprogram.equals("mm4")||qmprogram.equals("mm4hr")){
			qmMethod = "mm4";
			if(qmprogram.equals("mm4hr")) {
				qmMethod = "mm4hr";
				useHindRot=true;
			}
		}
		else{
			//may eventually want to pass this to various functions to choose which "sub-function" to call
			qmMethod="pm3";
		}

		mapMaxAttemptNumber = new HashMap<String, Integer>();
		mapMaxAttemptNumber.put("gaussian03", 36);
		mapMaxAttemptNumber.put("mopac", 10);
		mapMaxAttemptNumber.put("both", 10);
		mapMaxAttemptNumber.put("mm4", 4);
		mapMaxAttemptNumber.put("mm4hr", 4);


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
	/**
	 * //if there is no data in the libraries, calculate the result based on QM or MM calculations;
	 *  the below steps will be generalized later to allow for other quantum mechanics packages, etc.
	 *  
	 *  Several steps can be identified:
	 *  <LI> InChI generation of ChemGraph
	 *  <LI> Verification whether this species has already been processed (and parse it if so)
	 *  <LI> Generating 3D coords
	 *  <LI> Creating QM input file and Feeding species to QM Program
	 *  <LI> QM Program Output File Parsing
	 * @param p_chemGraph
	 * @return
	 */
	public ThermoData generateQMThermoData(ChemGraph p_chemGraph){

		//molfile that will store 3D coords
		molFile p_3dfile = null;

		ThermoData result = new ThermoData();
		double[] dihedralMinima = null;

		String [] InChInames = getQMFileName(p_chemGraph);//determine the filename (InChIKey) and InChI with appended info for triplets, etc.
		String name = InChInames[0];
		String InChIaug = InChInames[1];
		
		/*
		 * verify whether a succesful QM results exists for this particular species:
		 */
		QMVerifier verifier = new QMVerifier(name, InChIaug);
		verifier.verify();

		/*
		 * if a succesful job exists (by one of the QM Programs),
		 *  you can readily parse it.
		 */
		if(verifier.succesfulJobExists()){
			result = parseOutput(name, p_chemGraph);
			return result;
		}
		/*
		 * create 3D coords of ChemGraph:
		 */
		else{
			p_3dfile = generate3DCoords(p_chemGraph, name);


			int attemptNumber=1;//counter for attempts using different keywords

			//dynamic flag for success of Gaussian run; 0 means it failed, 1 means it succeeded
			int successFlag=0;

			/*
			 * get max number of attempts, based on the qm program that is used.
			 */
			int maxAttemptNumber = mapMaxAttemptNumber.get(qmprogram);
			int multiplicity = p_chemGraph.getRadicalNumber()+1; //multiplicity = radical number + 1

			while(successFlag==0 && attemptNumber <= maxAttemptNumber){

				createQMInput(name,p_chemGraph, p_3dfile, attemptNumber, InChIaug, multiplicity);

				successFlag = runQM(name);


				//new IF block to check success
				if(successFlag==1){
					Logger.info("Attempt #"+attemptNumber + " on species " + name + " ("+InChIaug+") succeeded.");
					//run rotor calculations if necessary
					int rotors = p_chemGraph.getInternalRotor();
					if(qmprogram.equals("mm4hr") && useHindRot && rotors > 0){
						//we should re-run scans even if pre-existing scans exist because atom numbering may change from case to case; a better solution would be to check for stored CanTherm output and use that if available
						Logger.info("Running rotor scans on "+name+"...");
						dihedralMinima = createMM4RotorInput(name, p_chemGraph, rotors);//we don't worry about checking InChI here; if there is an InChI mismatch it should be caught
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
							String directory = qmfolder;
							File dir=new File(directory);
							directory = dir.getAbsolutePath();//this and previous three lines get the absolute path for the directory
							QMJobRunnable job = new MM4HRJob(name, directory);
							job.run();
						}
					}
				}
				else if(successFlag==0){
					if(attemptNumber==maxAttemptNumber){//if this is the last possible attempt, and the calculation fails, exit with an error message
						if(qmprogram.equals("both")){ //if we are running with "both" option and all keywords fail, try with Gaussian
							qmprogram = "gaussian03";
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
			if(qmprogram.equals("mm4hr") && useCanTherm){
				performCanThermCalcs(name, p_chemGraph, dihedralMinima);
				if (p_chemGraph.getInternalRotor()>0) performCanThermCalcs(name, p_chemGraph, dihedralMinima);//calculate RRHO case for comparison
			}
			/*
			 * 5. parse QM output and record as thermo data (function includes symmetry/point group calcs,
			 *  etc.); if both Gaussian and MOPAC results exist, Gaussian result is used
			 */

			result = parseOutput(name, p_chemGraph);
			return result;
		}


	}
	private int runQM(String name) {
		String directory = qmfolder;
		File dir=new File(directory);
		directory = dir.getAbsolutePath();//this and previous three lines get the absolute path for the directory
		QMJobRunnable job;
		if (qmprogram.equals("gaussian03")){
			if(qmMethod.equals("pm3")){
				//4. run Gaussian
				/*
				 * name and directory are the name and directory for the input (and output) file;
				 * input is assumed to be preexisting and have the .gjf suffix
				 * returns an integer indicating success or failure of the Gaussian calculation: 1 for success, 0 for failure;
				 */
				job = new GaussianJob(name, directory);
				return job.run();
			}

		}

		else if(qmprogram.equals("mopac") || qmprogram.equals("both")){
			/*
			 * name and directory are the name and directory for the input (and output) file;
			 * input is assumed to be preexisting and have the .mop suffix
			 * returns an integer indicating success or failure of the MOPAC calculation: 1 for success, 0 for failure;
			 * this function is based on the Gaussian analogue
			 */
			job = new MOPACJob(name, directory); 
			return job.run();
		}
		else if(qmprogram.equals("mm4") || qmprogram.equals("mm4hr")){
			//4. run MM4
			//name and directory are the name and directory for the input (and output) file;
			//input script is assumed to be preexisting and have the .com suffix
			//returns an integer indicating success or failure of the calculation: 1 for success, 0 for failure
			job = new MM4Job(name, directory);
			return job.run();
		}
		else{
			Logger.critical("Unsupported quantum chemistry program");
			System.exit(0);
		}
		return -1;

	}
	private void createQMInput(String name, ChemGraph p_chemGraph, molFile p_3dfile, int attemptNumber,
			String inChIaug, int multiplicity) {
		//3. create the Gaussian or MOPAC input file
		String directory = qmfolder;
		File dir=new File(directory);
		directory = dir.getAbsolutePath();//this and previous three lines get the absolute path for the directory
		QMInputWritable writer;
		File inputFile;

		if (qmprogram.equals("gaussian03")){
			if(qmMethod.equals("pm3")){
				/**
				 *creates Gaussian PM3 input file in directory with filename name.gjf by using OpenBabel
				 * to convert p_molfile
				 *@param attemptNumber determines which keywords to try
				 *@return the function returns the maximum number of keywords that can be attempted; 
				 *this will be the same throughout the evaluation of the code, 
				 *so it may be more appropriate to have this as a "constant" attribute of some sort
				 *attemptNumber=-1 will call a special set of keywords for the monoatomic case     
				 */
				if(p_chemGraph.getAtomNumber() > 1){
					//write a file with the input keywords
					writer = new GaussianPM3InputWriter(name, directory, p_3dfile, attemptNumber, multiplicity, inChIaug);

					inputFile = writer.write();
				}
				//use -1 for attemptNumber for monoatomic case
				else{
					//write a file with the input keywords
					writer = new GaussianPM3InputWriter(name, directory, p_3dfile, -1, multiplicity, inChIaug);
					inputFile = writer.write();

				}
			}
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
		else if(qmprogram.equals("mopac") || qmprogram.equals("both")){
			//write a file with the input keywords
			writer = new MOPACPM3InputWriter(name, directory, p_3dfile, attemptNumber, multiplicity);
			inputFile = writer.write();
		}
		else if(qmprogram.equals("mm4") || qmprogram.equals("mm4hr")){
			//write a file with the input keywords
			writer = new MM4InputWriter(name, directory, p_3dfile, attemptNumber, multiplicity, inChIaug);

			inputFile = writer.write();

		}


	}
	public molFile generate3DCoords(ChemGraph p_chemGraph, String name) {
		molFile p_3dfile;
		//steps 1 and 2: create 2D and 3D mole files
		ThreeDMolFileCreator creator = new ThreeDMolFileCreator(name, p_chemGraph);
		p_3dfile = creator.create();
		return p_3dfile;
	}

	/**
	 * wrapper method for parser types
	 * @param name
	 * @param p_chemGraph
	 * @return
	 */
	public ThermoData parseOutput(String name, ChemGraph p_chemGraph) {
		String directory = qmfolder;
		File dir=new File(directory);
		directory = dir.getAbsolutePath();//this and previous three lines get the absolute path for the directory

		ThermoData result = null;
		if(qmMethod.equals("pm3")){
			if ((qmprogram.equals("gaussian03"))){
				//parse the results using cclib and return a ThermoData object; name and directory indicate the location of the Gaussian .log file
				//may want to split this into several functions
				QMParsable parser = new GaussianPM3Parser(name, directory, p_chemGraph);
				result = parser.parse();
				result.setSource("Gaussian PM3 calculation");
				Logger.info("Thermo for " + name + ": "+ result.toString());//print result, at least for debugging purposes
				return result;

			}
			else if (qmprogram.equals("mopac") || qmprogram.equals("both")){
				QMParsable parser = new MOPACPM3Parser(name, directory, p_chemGraph);

				result = parser.parse();
				result.setSource("MOPAC PM3 calculation");
				Logger.info("Thermo for " + name + ": "+ result.toString());//print result, at least for debugging purposes
				return result;
			}
			else if (qmprogram.equals("mm4")){
				//5. parse MM4 output and record as thermo data (function includes symmetry/point group calcs, etc.)
				if(!useCanTherm) {
					QMParsable parser = new MM4Parser(name, directory, p_chemGraph);

					result = parser.parse();

					result.setSource("MM4 calculation");
					Logger.info("Thermo for " + name + ": "+ result.toString());//print result, at least for debugging purposes
				}
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
			else{
				Logger.critical("Unexpected situation in QMTP thermo estimation");
				System.exit(0);
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
	 * creates MM4 rotor input file and MM4 batch file in directory 
	 * with filenames name.mm4roti and name.comi, respectively
	 * the function returns the set of rotor dihedral angles for 
	 * the minimum energy conformation 
	 */
	public double[] createMM4RotorInput(String name, ChemGraph p_chemgraph, int rotors){
		String directory = qmfolder;
		File dir=new File(directory);
		directory = dir.getAbsolutePath();//this and previous three lines get the absolute path for the directory

		
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

	//parse the results using cclib and CanTherm and return a ThermoData object; name and directory indicate the location of the MM4 .mm4out file
	//formerly known as parseMM4withForceMat
	public IQMData performCanThermCalcs(String name, ChemGraph p_chemGraph, double[] dihedralMinima){
		String directory = qmfolder;
		File dir=new File(directory);
		directory = dir.getAbsolutePath();//this and previous three lines get the absolute path for the directory

		
		// parse the MM4 file with cclib to get atomic number vector and geometry
		QMParser parser = new MM4Parser(name, directory, p_chemGraph, true);

		/*
		 * TODO if we don't call the parse() method the IQMData object
		 * would not be created...
		 */
		ThermoData data = parser.parse();

		IQMData qmdata = parser.getQMData();
		//unpack the needed results and write cantherm input file
		//flag RRHO = true
		QMInputWritable canthermWriter = new CanThermInputWriter(name, directory, p_chemGraph, qmdata, dihedralMinima, true);
		File canThermFile = canthermWriter.write();

		// call CanTherm 
		QMJobRunnable canthermJob = new CanThermJob(name, directory);
		canthermJob.run();

		return qmdata;
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
