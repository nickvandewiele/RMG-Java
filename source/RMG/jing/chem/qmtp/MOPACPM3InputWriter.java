package jing.chem.qmtp;

import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Map;

import jing.chem.QMTP;
import jing.chem.molFile;
import jing.rxnSys.Logger;

public class MOPACPM3InputWriter extends QMInputWriter implements QMInputWritable{

	public Map<Integer, String> multiplicityKeywords = new HashMap<Integer, String>();
	
	public String InChIAug;
	
	private class MOPACKEYWORDS{
		public static final String BOTH = "Both";
		public static final String BOTTOM = "Bottom";
		public static final String TOP = "Top";
		public static final String POLAR = "Polar";
	}
	public int multiplicity;
	public MOPACPM3InputWriter(
			String name, 
			String directory,
			molFile p_molfile, 
			int attemptNumber, 
			int multiplicity
			) {
		
		super(name, directory, p_molfile, attemptNumber, multiplicity);
		
		scriptAttempts = 5;
		
		maxAttemptNumber = 2 * scriptAttempts;
		
		fillMultiplicityKeywords();
	}

	private void fillMultiplicityKeywords() {
		
		multiplicityKeywords.put(1, "");
		multiplicityKeywords.put(2, "uhf doublet");
		multiplicityKeywords.put(3, "uhf triplet");
		multiplicityKeywords.put(4, "uhf quartet");
		multiplicityKeywords.put(5, "uhf quintet");
		multiplicityKeywords.put(6, "uhf sextet");
		multiplicityKeywords.put(7, "uhf septet");
		multiplicityKeywords.put(8, "uhf octet");
		multiplicityKeywords.put(9, "uhf nonet");
		
	}

	public File write() {
		createKeywords();

		File inputFile = createInputFile();
		
		return inputFile;
	}
	/**
	 * returns the extra Mopac keywords to use for radical species, given the spin multiplicity 
	 * (radical number + 1)
	 * @param multipl
	 * @return
	 */
	public String getMopacRadicalString(int multipl){
		
		String keyword = null;
		try{
			keyword = multiplicityKeywords.get(new Integer(multipl));	
		}
		catch(Exception e){
			Logger.critical("Invalid multiplicity encountered: "+multipl);
			Logger.critical("this should not be returned: error associated with getMopacRadicalString()");
			System.exit(0);
		}

		return keyword;

	}
	/**
	 *  call the OpenBabel process (note that this requires OpenBabel environment variable)
	 * @return
	 */
	public File createInputFile() {
		//call the OpenBabel process (note that this requires OpenBabel environment variable)
		try{ 
			File runningdir=new File(dir);
			String inpKeyStrTopCombined = keywords.get("Both") + keywords.get("Top");
			Process babelProc = null;
			if(attemptNumber<=scriptAttempts){//use UFF-refined coordinates
				String[] command = { "babel", "-imol", molfile.getPath(), "-xk", inpKeyStrTopCombined,"--title", InChIAug,"-omop", name+".mop" };
				babelProc = Runtime.getRuntime().exec(command, null, runningdir);
			}
			else{//use the crude coordinates
				String[] command = { "babel", "-imol", molfile.getCrudePath(), "-xk",  inpKeyStrTopCombined,"--title", InChIAug,"-omop", name+".mop" };
				babelProc = Runtime.getRuntime().exec(command, null, runningdir);
			}


			//append the final keywords to the end of the file just written
			// File mopacInpFile = new File(directory+"/"+name+".mop");
			FileWriter fw = new FileWriter(dir+"/"+name+".mop", true);//filewriter with append = true
			fw.write(System.getProperty("line.separator") + keywords.get("Bottom") + keywords.get("Both") + keywords.get("Polar"));//on Windows Vista, "\n" appears correctly in WordPad view, but not Notepad view (cf. http://forums.sun.com/thread.jspa?threadID=5386822)
			fw.close();
		}
		catch(Exception e){
			String err = "Error in running OpenBabel MOL to MOP process \n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}
		return new File(name+".mop");
	}

	@Override
	public Map<String, String> createKeywords() {
		/**
		 * this string will be written at both the top (for optimization) and the bottom (for thermo/force calc)
		 */
		String inpKeyStrBoth = "";
		/**
		 * this string will be written only at the top
		 */
		String inpKeyStrTop = "";
		/**
		 * this string will be written at the bottom
		 */
		String inpKeyStrBottom = "";
		
		String radicalString = getMopacRadicalString(multiplicity);
		try{
			//      File inpKey=new File(directory+"/inputkeywords.txt");
			//      String inpKeyStr="%chk="+directory+"\\RMGrunCHKfile.chk\n";
			//      inpKeyStr+="%mem=6MW\n";
			//     inpKeyStr+="%nproc=1\n";
			switch(attemptNumber%scriptAttempts){
			case 1:
				inpKeyStrBoth="pm3 "+radicalString;
				inpKeyStrTop=" precise nosym";
				/**
				 * 7/10/09: based on a quick review of recent results, 
				 * keyword combo #1 rarely works, and when it did (CJAINEUZFLXGFA-UHFFFAOYAUmult3 
				 * (InChI=1/C8H16O5Si/c1-4-11-14(9,12-5-2)13-8-6-10-7(8)3/h7-8H,3-6H2,1-2H3/mult3)), 
				 * the grad. norm on the force step was about 1.7 (too large); 
				 * I manually removed this result and re-ran...the entropy was increased 
				 * by nearly 20 cal/mol-K...perhaps we should add a check for the "WARNING" 
				 * that MOPAC prints out when the gradient is high; 7/22/09: for the case of 
				 * FUGDBSHZYPTWLG-UHFFFAOYADmult3 (InChI=1/C5H8/c1-4-3-5(4)2/h4-5H,1-3H2/mult3), 
				 * adding nosym seemed to resolve 1. large discrepancies from Gaussian and 2. 
				 * negative frequencies in mass-weighted coordinates and possibly related issue 
				 * in discrepancies between regular and mass-weighted coordinate frequencies
				 */
				inpKeyStrBottom="oldgeo thermo nosym precise ";
				break;
			case 2:
				/**
				 * //7/9/09: used for VCSJVABXVCFDRA-UHFFFAOYAI 
				 * (InChI=1/C8H19O5Si/c1-5-10-8(4)13-14(9,11-6-2)12-7-3/h8H,5-7H2,1-4H3); 
				 * all existing Gaussian keywords also failed; the Gaussian result 
				 * was also rectified, but the resulting molecule was over 70 kcal/mol less stable,
				 *  probably due to a large amount of spin contamination (~1.75 in fixed Gaussian 
				 *  result vs. 0.754 for MOPAC)
				 */
				inpKeyStrBoth="pm3 "+radicalString;
				inpKeyStrTop=" precise nosym gnorm=0.0 nonr";
				inpKeyStrBottom="oldgeo thermo nosym precise ";
				break;
				
			case 3:
				/**
				 * //7/8/09: used for ADMPQLGIEMRGAT-UHFFFAOYAUmult3 
				 * (InChI=1/C6H14O5Si/c1-4-9-12(8,10-5-2)11-6(3)7/h6-7H,3-5H2,1-2H3/mult3); 
				 * all existing Gaussian keywords also failed; however, the Gaussian result
				 *  was also rectified, and the resulting conformation was about 1.0 kcal/mol 
				 *  more stable than the one resulting from this, so fixed Gaussian result was 
				 *  manually copied to QMFiles folder
				 */
				inpKeyStrBoth="pm3 "+radicalString;
				inpKeyStrTop=" precise nosym gnorm=0.0";
				/**
				 * //precise appeared to be necessary for the problematic case (to avoid negative frequencies);
				 */
				inpKeyStrBottom="oldgeo thermo nosym precise "; 
				break;
			case 4:
				/**
				 * //7/8/09: used for GYFVJYRUZAKGFA-UHFFFAOYALmult3 
				 * (InChI=1/C6H14O6Si/c1-3-10-13(8,11-4-2)12-6-5-9-7/h6-7H,3-5H2,1-2H3/mult3) 
				 * case (negative frequency issues in MOPAC) (also, none of the existing Gaussian
				 *  combinations worked with it); note that the Gaussian result appears to be a 
				 *  different conformation as it is about 0.85 kcal/mol more stable, so the Gaussian 
				 *  result was manually copied to QMFiles directory; note that the MOPAC output included 
				 *  a very low frequency (4-5 cm^-1)
				 */
				inpKeyStrBoth="pm3 "+radicalString;
				inpKeyStrTop=" precise nosym gnorm=0.0 bfgs";
				/**
				 * precise appeared to be necessary for the problematic case (to avoid negative frequencies)
				 */
				inpKeyStrBottom="oldgeo thermo nosym precise "; 
				break;
				//            else if(attemptNumber==5){
				//                inpKeyStrBoth="pm3 "+radicalString;
				//                inpKeyStrTop=" precise nosym gnorm=0.0 ddmin=0.0";
				//                inpKeyStrBottom="oldgeo thermo nosym precise ";
				//            }
				//            else if(attemptNumber==6){
				//                inpKeyStrBoth="pm3 "+radicalString;
				//                inpKeyStrTop=" precise nosym gnorm=0.0 nonr ddmin=0.0";
				//                inpKeyStrBottom="oldgeo thermo nosym precise ";
				//            }
				//            else if(attemptNumber==7){
				//                inpKeyStrBoth="pm3 "+radicalString;
				//                inpKeyStrTop=" precise nosym bfgs gnorm=0.0 ddmin=0.0";
				//                inpKeyStrBottom="oldgeo thermo nosym precise ";
				//            }
			case 0:
				/**
				 * used for troublesome HGRZRPHFLAXXBT-UHFFFAOYAVmult3 
				 * (InChI=1/C3H2O4/c4-2(5)1-3(6)7/h1H2/mult3) case 
				 * (negative frequency and large gradient issues)				
				 */
				inpKeyStrBoth="pm3 "+radicalString;
				inpKeyStrTop=" precise nosym recalc=10 dmax=0.10 nonr cycles=2000 t=2000";
				inpKeyStrBottom="oldgeo thermo nosym precise ";
				//   else if(attemptNumber==9){//used for troublesome CMARQPBQDRXBTN-UHFFFAOYAAmult3 (InChI=1/C3H2O4/c1-3(5)7-6-2-4/h1H2/mult3) case (negative frequency issues)
				//       inpKeyStrBoth="pm3 "+radicalString;
				//       inpKeyStrTop=" precise nosym recalc=1 dmax=0.05 gnorm=0.0 cycles=1000 t=1000";
				//       inpKeyStrBottom="oldgeo thermo nosym precise ";
				//   }
				//   else if(attemptNumber==10){//used for ATCYLHQLTOSVFK-UHFFFAOYAMmult4 (InChI=1/C4H5O5/c1-3(5)8-9-4(2)7-6/h6H,1-2H2/mult4) case (timeout issue; also, negative frequency issues); note that this is very similar to the keyword below, so we may want to consolidate
				//       inpKeyStrBoth="pm3 "+radicalString;
				//       inpKeyStrTop=" precise nosym recalc=1 dmax=0.05 gnorm=0.2 cycles=1000 t=1000";
				//       inpKeyStrBottom="oldgeo thermo nosym precise ";
				//   }
				break;
			default:
				throw new Exception();//this point should not be reached
				
			}

		}
		catch(Exception e){
			String err = "Error in writing inputkeywords.txt \n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}

		String polarString = "";
		if (QMTP.usePolar){
			if(multiplicity == 1) polarString = System.getProperty("line.separator") + System.getProperty("line.separator") + System.getProperty("line.separator")+ "oldgeo polar nosym precise " + inpKeyStrBoth;
			else polarString = System.getProperty("line.separator") + System.getProperty("line.separator") + System.getProperty("line.separator")+ "oldgeo static nosym precise " + inpKeyStrBoth;
		}
		
		keywords = new HashMap<String, String>();
		keywords.put(MOPACKEYWORDS.BOTH, inpKeyStrBoth);
		keywords.put(MOPACKEYWORDS.TOP, inpKeyStrTop);
		keywords.put(MOPACKEYWORDS.BOTTOM, inpKeyStrBottom);
		keywords.put(MOPACKEYWORDS.POLAR, polarString);
		
		return keywords;
	}



}
