package jing.chem.qmtp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;

import jing.chem.ChemGraph;
import jing.chem.QMTP;
import jing.chem.molFile;
import jing.rxnSys.Logger;

public class MM4HRInputWriter extends QMInputWriter implements QMInputWritable {

	public String InChIAug;
	
	public ChemGraph chemGraph;
	
	public int indexForFirstAtom;
	
	String[] lines;
	
	String mm4optContents;
	
	double[] dihedralMinima;
	
	private class MM4HRKEYWORDS{
		public static final String INPUT = "Input";
	}
	
	public MM4HRInputWriter(
			String name, 
			String directory, 
			ChemGraph chemGraph,
			int rotors

	) {

		super(name, directory);
		this.chemGraph = chemGraph;
		scriptAttempts = 2;

		maxAttemptNumber = 2 * scriptAttempts;

		indexForFirstAtom=-1;
		
		mm4optContents = null;
		
		dihedralMinima = new double[rotors];
		
	}

	@Override
	public File write() {
		createInputFile();
		/*
		 * for now, null will be returned,
		 * returning one single file does not make sense since multiple files are created here...
		 */
		return null;
	}

	@Override
	/**
	 * this method is not used because each single file containing a specific rotation of
	 * the rotor has a specific index i attached to it.
	 * 
	 * So, instead I placed the keyword generation procedure directly in the 
	 * {@link MM4HRInputWriter.createInputFile()} method.
	 */
	public Map<String, String> createKeywords() {
		
		return null;

	}

	@Override
	public File createInputFile() {
		//iterate over all the rotors in the molecule
		int i = 0;//rotor index
		Map rotorInfo = chemGraph.getInternalRotorInformation();
		Iterator iter = rotorInfo.keySet().iterator();
		while(iter.hasNext()){
			i++;
			int[] rotorAtoms = (int[])iter.next();
			String mm4optContents = readMM4OptFileContents();
			String inpKeyStr = null;
			/*
			 * split by newlines, excluding blank lines; 
			 * cf. http://stackoverflow.com/questions/454908/split-java-string-by-new-line
			 */
			String[] lines = mm4optContents.split("[\\r?\\n]+");
			//flag to indicate whether there is a "line 1a" with pi system information
			int flag = 0;
			//pi system information is indicated by the second line beginning with T's and potentially F's
			if (lines[1].startsWith("T")||lines[1].startsWith("F")) flag = 1;
			/*
			 * take the first 78 characters of line 2 (or 3 if it is a pi-system compound), 
			 * and append the option number for the NDRIVE option; in other words, 
			 * we are replacing the NDRIVE=0 option with the desired option number
			 */
			lines[1+flag]=lines[1+flag].substring(0,78)+" 2";
			//reconstruct mm4optContents
			mm4optContents = "";
			for(int j=0; j<lines.length;j++){
				mm4optContents +=lines[j]+"\n";
			}


			try{
				//write one script file for each rotor
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
				/*
				 * create batch file with executable permissions:
				 * cf. http://java.sun.com/docs/books/tutorial/essential/io/fileAttr.html#posix
				 */

				inpKeyStr="#! /bin/csh\n";
				inpKeyStr+="cp "+name+".mm4rot"+i+" CPD.MM4\n";
				inpKeyStr+="cp $MM4_DATDIR/BLANK.DAT PARA.MM4\n";
				inpKeyStr+="cp $MM4_DATDIR/CONST.MM4 .\n";
				inpKeyStr+="$MM4_EXEDIR/mm4 <<%\n";
				inpKeyStr+="1\n";//read from first line of .mm4 file
				inpKeyStr+="3\n"; //Full-Matrix Method only
				//inpKeyStr+="2\n"; //Block-Diagonal Method then Full-Matrix Method

				inpKeyStr+="0\n";//terminate the job
				inpKeyStr+="%\n";
				inpKeyStr+="mv TAPE4.MM4 "+name+".mm4rotout"+i+"\n";
				inpKeyStr+="mv TAPE9.MM4 "+name+".mm4rotopt"+i+"\n";
				inpKeyStr+="exit\n";

			}
			catch(Exception e){
				String err = "Error in writing MM4 script files\n";
				err += e.toString();
				Logger.logStackTrace(e);
				System.exit(0);
			}

			/*
			 * keywords attribute is always re-initialized depending on the specific 
			 * rotation of the rotor at hand...µ
			 * thus, keywords should not be an attribute in the long run.
			 */
			keywords = new HashMap<String, String>();
			keywords.put(MM4HRKEYWORDS.INPUT, inpKeyStr);

			File inpKey = new File(dir+"/"+name+".com"+i);
			FileWriter fw;
			try {
				fw = new FileWriter(inpKey);
				fw.write(keywords.get(MM4HRKEYWORDS.INPUT));
				fw.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			//Step 2: create the MM4 rotor job input file for rotor i, using the output file from the "normal" MM4 job as a template
			//Step 2a: find the dihedral angle for the minimized conformation (we need to start at the minimum energy conformation because CanTherm assumes that theta=0 is the minimum
			//extract the lines for dihedral atoms
			readLines();
			String dihedral1s = lines[indexForFirstAtom+rotorAtoms[0]-1];
			String atom1s = lines[indexForFirstAtom+rotorAtoms[1]-1];
			String atom2s = lines[indexForFirstAtom+rotorAtoms[2]-1];
			String dihedral2s = lines[indexForFirstAtom+rotorAtoms[3]-1];
			//extract the x,y,z coordinates for each atom
			double[] dihedral1 = {
					Double.parseDouble(dihedral1s.substring(0,10)),
					Double.parseDouble(dihedral1s.substring(10,20)),
					Double.parseDouble(dihedral1s.substring(20,30))};
			
			double[] atom1 = {
					Double.parseDouble(atom1s.substring(0,10)),
					Double.parseDouble(atom1s.substring(10,20)),
					Double.parseDouble(atom1s.substring(20,30))
					};
			
			double[] atom2 = {
					Double.parseDouble(atom2s.substring(0,10)),
					Double.parseDouble(atom2s.substring(10,20)),
					Double.parseDouble(atom2s.substring(20,30))
					};
			
			double[] dihedral2 = {
					Double.parseDouble(dihedral2s.substring(0,10)),
					Double.parseDouble(dihedral2s.substring(10,20)),
					Double.parseDouble(dihedral2s.substring(20,30))
					};
			//determine the dihedral angle
			dihedralMinima[i-1] = calculateDihedral(dihedral1,atom1,atom2,dihedral2);
			
			//double dihedral = calculateDihedral(dihedral1,atom1,atom2,dihedral2);
			//if (dihedral < 0) dihedral = dihedral + 360;
			/*
			 * make sure dihedral is positive; this will save 
			 * an extra character in the limited space for specifying starting and ending angles
			 * 
			 * eventually, when problems arise due to collinear atoms 
			 * (both arguments to atan2 are zero) we can iterate over 
			 * other atom combinations (with atoms in each piece determined 
			 * by the corresponding value in rotorInfo for atom2 and the complement
			 *  of these (full set = p_chemgraph.getNodeIDs()) for atom1) until they are not collinear
			 */

			//Step 2b: write the file for rotor i
			try{
				FileWriter mm4roti = new FileWriter(dir+"/"+name+".mm4rot"+i);
				/*
				 * //deltaTheta should be less than 100 degrees so that dihedral-deltaTheta still fits
				 */
				//mm4roti.write(mm4optContents+"\n"+String.format("  %3d  %3d  %3d  %3d     %5.1f%5.1f%5.1f", rotorAtoms[0],rotorAtoms[1],rotorAtoms[2],rotorAtoms[3], dihedral, dihedral + 360.0 - deltaTheta, deltaTheta)+"\n");
				//deltaTheta should be less than 100 degrees so that dihedral-deltaTheta still fits
				//mm4roti.write(mm4optContents+"\n"+String.format("  %3d  %3d  %3d  %3d     %5.1f%5.1f%5.1f", rotorAtoms[0],rotorAtoms[1],rotorAtoms[2],rotorAtoms[3], dihedral, dihedral - deltaTheta, deltaTheta)+"\n");
				/*
				 * //M1, M2, M3, M4, START, FINISH, DIFF (as described in MM4 manual); 
				 * //this would require miminum to be stored and used to adjust 
				 * actual angles before sending to CanTherm
				 */
				mm4roti.write(mm4optContents+String.format("  %3d  %3d  %3d  %3d     %5.1f%5.1f%5.1f", rotorAtoms[0],rotorAtoms[1],rotorAtoms[2],rotorAtoms[3], 0.0, 360.0-QMTP.deltaTheta, QMTP.deltaTheta)+"\n");
				mm4roti.close();
			}
			catch(Exception e){
				String err = "Error in writing MM4 rotor input file\n";
				err += e.toString();
				Logger.logStackTrace(e);
				System.exit(0);
			}
		}
		/*
		 * for now, null will be returned,
		 * returning one single file does not make sense since multiple files are created here...
		 */
		return null;
	}
	/**
	 * given x,y,z (cartesian) coordinates for dihedral1 atom, (rotor) atom 1, 
	 * (rotor) atom 2, and dihedral2 atom, calculates the dihedral angle 
	 * (in degrees, between +180 and -180) using the atan2 formula at 
	 * http://en.wikipedia.org/w/index.php?title=Dihedral_angle&oldid=373614697
	 * @param dihedral1
	 * @param atom1
	 * @param atom2
	 * @param dihedral2
	 * @return
	 */
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
	 * read in the optimized coordinates from the completed "normal" MM4 job; 
	 * this will be used as a template for the rotor input files
	 * @return
	 */
	public String readMM4OptFileContents(){

		mm4optContents = "";
		
		int lineIndex=0;
		try{
			FileReader mm4opt = new FileReader(dir+"/"+name+".mm4opt");
			BufferedReader reader = new BufferedReader(mm4opt);
			String line=reader.readLine();
			while(line!=null){
				if(indexForFirstAtom < 0 && line.length()>=40 && line.substring(35,40).equals("(  1)")) {
					indexForFirstAtom = lineIndex;//store the location of the first atom coordinate
				}
				mm4optContents += line + "\n";
				line=reader.readLine();
				lineIndex++;
			}
			reader.close();
			mm4opt.close();
		}
		catch(Exception e){
			String err = "Error in reading .mm4opt file\n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}
		return mm4optContents;
	}
	public void readLines(){
		/*
		 * split by newlines, excluding blank lines; 
		 * cf. http://stackoverflow.com/questions/454908/split-java-string-by-new-line
		 */
		lines = mm4optContents.split("[\\r?\\n]+");
		
		//flag to indicate whether there is a "line 1a" with pi system information
		int flag = 0;
		
		//pi system information is indicated by the second line beginning with T's and potentially F's
		if (lines[1].startsWith("T")||lines[1].startsWith("F")) flag = 1;
		
		/*
		 *  take the first 78 characters of line 2
		 *  (or 3 if it is a pi-system compound), 
		 *  and append the option number for the NDRIVE option; 
		 *  in other words, we are replacing the NDRIVE=0 option with the desired option number
		 */
		lines[1+flag]=lines[1+flag].substring(0,78)+" 2";
	}

	public double[] getDihedralMinima() {
		return dihedralMinima;
	}
}
