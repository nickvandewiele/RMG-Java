package jing.chem.qmtp;

import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Map;

import jing.chem.QMTP;
import jing.chem.molFile;
import jing.rxnSys.Logger;

public class GaussianPM3InputWriter extends QMInputWriter implements QMInputWritable{

	
	
	public String InChIAug;
	
	private class G03PM3KEYWORDS{
		public static final String INPUT = "Input";
	}

	public GaussianPM3InputWriter(String name, String directory,
			molFile p_molfile, int attemptNumber, int multiplicity,
			String inChIaug) {
		
		super(name, directory, p_molfile, attemptNumber, multiplicity);
		
		scriptAttempts = 18;
		
		maxAttemptNumber = 2 * scriptAttempts;
		
		
		this.InChIAug = inChIaug;
		
	}

	public File write() {
		createKeywords();

		File inputFile = createInputFile();
		
		return inputFile;
	}

	/**
	 *  call the OpenBabel process (note that this requires OpenBabel environment variable)
	 * @return
	 */
	public File createInputFile() {
		try{ 

			String command=null;
			String molPath=null;
			if(attemptNumber<= scriptAttempts){//use UFF-refined coordinates
				molPath = molfile.getPath();
			}
			else{//use crude coordinates
				molPath = molfile.getCrudePath();
			}
			if (System.getProperty("os.name").toLowerCase().contains("windows")){//special windows case
				command = "babel -imol \""+ molPath+ "\" -ogjf \"" + name+".gjf\" -xf inputkeywords.txt --title \""+InChIAug+"\"";
			}
			else{
				command = "babel -imol "+ molPath+ " -ogjf " + name+".gjf -xf inputkeywords.txt --title "+InChIAug;
			}
			
			Process babelProc = Runtime.getRuntime().exec(command, null, new File(dir));

			
		}
		catch(Exception e){
			String err = "Error in running OpenBabel MOL to GJF process \n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}
		

		return new File(name + ".gjf");
	}


	@Override
	public Map<String, String> createKeywords() {
		String inpKeyStr = null;
		try{
			File inputKeywords=new File(dir+"/inputkeywords.txt");
			inpKeyStr="%chk="+dir+"/RMGrunCHKfile.chk\n";
			inpKeyStr+="%mem=6MW\n";
			inpKeyStr+="%nproc=1\n";
			if(attemptNumber==-1) inpKeyStr+="# pm3 freq";//keywords for monoatomic case (still need to use freq keyword to get molecular mass)

			switch(attemptNumber%scriptAttempts){
			case 1:
				/**
				 * added IOP option to avoid aborting when symmetry changes; 
				 * 3 is supposed to be default according to documentation, 
				 * but it seems that 0 (the default) is the only option that doesn't work from 0-4;
				 *  also, it is interesting to note that all 4 options seem to work for test case 
				 *  with z-matrix input rather than xyz coords; 
				 *  cf. http://www.ccl.net/cgi-bin/ccl/message-new?2006+10+17+005 for original idea for 
				 *  solution
				 */
				inpKeyStr+="# pm3 opt=(verytight,gdiis) freq IOP(2/16=3)";
				break;
			case 2:
				/**
				 * use different SCF method; this addresses at least one case of failure for a C4H7J species
				 */
				inpKeyStr+="# pm3 opt=(verytight,gdiis) freq IOP(2/16=3) IOP(4/21=2)";
				break;

			case 3:
				/**
				 * try multiple different options (no gdiis, use calcfc, nosymm); 
				 * 7/21/09: added maxcyc option to fix case of MPTBUKVAJYJXDE-UHFFFAOYAPmult3 
				 * (InChI=1/C4H10O5Si/c1-3-7-9-10(5,6)8-4-2/h4-5H,3H2,1-2H3/mult3) 
				 * (file manually copied to speed things along)
				 */

				inpKeyStr+="# pm3 opt=(verytight,calcfc,maxcyc=200) freq IOP(2/16=3) nosymm";
				break;
			case 4:
				/**
				 * //7/8/09: numerical frequency keyword version of keyword #3; 
				 * used to address GYFVJYRUZAKGFA-UHFFFAOYALmult3 
				 * (InChI=1/C6H14O6Si/c1-3-10-13(8,11-4-2)12-6-5-9-7/h6-7H,3-5H2,1-2H3/mult3) case; 
				 * (none of the existing Gaussian or MOPAC combinations worked with it)
				 */
				inpKeyStr+="# pm3 opt=(verytight,calcfc,maxcyc=200) freq=numerical IOP(2/16=3) nosymm";
				break;
			case 5:
				/**
				 * 7/10/09: somehow, this worked for problematic case of ZGAWAHRALACNPM-UHFFFAOYAF 
				 * (InChI=1/C8H17O5Si/c1-3-11-14(10,12-4-2)13-8-5-7(9)6-8/h7-9H,3-6H2,1-2H3); 
				 * (was otherwise giving l402 errors); even though I had a keyword that worked for this case,
				 *  I manually copied the fixed log file to QMfiles folder to speed things along;
				 *   note that there are a couple of very low frequencies (~5-6 cm^-1 for this case)
				 */
				inpKeyStr+="# pm3 opt=(verytight,gdiis,small) freq IOP(2/16=3)";
				break;
			case 6:
				/**
				 * //used for troublesome C5H7J2 case (similar error to C5H7J below); 
				 * calcfc is not necessary for this particular species,
				 *  but it speeds convergence and probably makes it more robust for other species
				 */
				inpKeyStr+="# pm3 opt=(verytight,nolinear,calcfc,small) freq IOP(2/16=3)";
				break;
			case 7:
				/**
				 * use numerical frequencies; this takes a relatively long time, 
				 * so should only be used as one of the last resorts; 
				 * this seemed to address at least one case of failure for a C6H10JJ species; 
				 * 7/15/09: maxcyc=200 added to address GVCMURUDAUQXEY-UHFFFAOYAVmult3 
				 * (InChI=1/C3H4O7Si/c1-2(9-6)10-11(7,8)3(4)5/h6-7H,1H2/mult3)...however, 
				 * result was manually pasted in QMfiles folder to speed things along
				 */
				inpKeyStr+="# pm3 opt=(verytight,gdiis,maxcyc=200) freq=numerical IOP(2/16=3)"; 
				break;
			case 8:
				/**
				 * //7/10/09: this worked for problematic case of SZSSHFMXPBKYPR-UHFFFAOYAF 
				 * (InChI=1/C7H15O5Si/c1-3-10-13(8,11-4-2)12-7-5-6-9-7/h7H,3-6H2,1-2H3) 
				 * (otherwise, it had l402.exe errors); corrected log file was manually copied to 
				 * QMfiles to speed things along; we could also add a freq=numerical 
				 * version of this keyword combination for added robustness; UPDATE: see below
				 */
				inpKeyStr+="# pm3 opt=tight freq IOP(2/16=3)";
				break;
			case 9:
				/**
				 * 7/10/09: used for problematic case of CIKDVMUGTARZCK-UHFFFAOYAImult4 
				 * (InChI=1/C8H15O6Si/c1-4-12-15(10,13-5-2)14-7-6-11-8(7,3)9/h7H,3-6H2,1-2H3/mult4 
				 * (most other cases had l402.exe errors);  corrected log file was manually copied to 
				 * QMfiles to speed things along
				 */
				inpKeyStr+="# pm3 opt=tight freq=numerical IOP(2/16=3)";
				break;
			case 10:
				/**
				 * 7/8/09: similar to existing #5, but uses tight rather than verytight; 
				 * used for ADMPQLGIEMRGAT-UHFFFAOYAUmult3 
				 * (InChI=1/C6H14O5Si/c1-4-9-12(8,10-5-2)11-6(3)7/h6-7H,3-5H2,1-2H3/mult3)
				 */
				inpKeyStr+="# pm3 opt=(tight,nolinear,calcfc,small,maxcyc=200) freq IOP(2/16=3)";
				break;
			case 11:
				/**
				 * use default (not verytight) convergence criteria; use this as last resort
				 */
				inpKeyStr+="# pm3 opt freq IOP(2/16=3)"; 
				break;
			case 12:
				/**
				 * to address problematic C10H14JJ case
				 */
				inpKeyStr+="# pm3 opt=(verytight,gdiis) freq=numerical IOP(2/16=3) IOP(4/21=200)";
				break;
			case 13:
				/**
				 * added 6/10/09 for very troublesome RRMZRNPRCUANER-UHFFFAOYAQ
				 *  (InChI=1/C5H7/c1-3-5-4-2/h3H,1-2H3) case...there were troubles with negative frequencies,
				 *   where I don't think they should have been; step size of numerical frequency was 
				 *   adjusted to give positive result; accuracy of result is questionable; 
				 *   it is possible that not all of these keywords are needed; note that for this
				 *    and other nearly free rotor cases, I think heat capacity will be 
				 *    overestimated by R/2 (R vs. R/2) (but this is a separate issue)
				 */
				inpKeyStr+="# pm3 opt=(calcfc,verytight,newton,notrustupdate,small,maxcyc=100,maxstep=100) freq=(numerical,step=10) IOP(2/16=3) nosymm";
				break;
			case 14:
				/**
				 * added 6/22/09 for troublesome QDERTVAGQZYPHT-UHFFFAOYAHmult3
				 * (InChI=1/C6H14O4Si/c1-4-8-11(7,9-5-2)10-6-3/h4H,5-6H2,1-3H3/mult3);
				 * key aspects appear to be tight (rather than verytight) convergence criteria, 
				 * no calculation of frequencies during optimization, use of numerical frequencies,
				 * and probably also the use of opt=small
				 */

				inpKeyStr+="# pm3 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm";
				break;
			case 15:
				/**
				 * gmagoon 7/9/09: commented out since although this produces a 
				 * "reasonable" result for the problematic case, there is a large amount 
				 * of spin contamination, apparently introducing 70+ kcal/mol of instability    
				 * else if(attemptNumber==12) inpKeyStr+="# pm3 opt=(verytight,gdiis,small) 
				 * freq=numerical IOP(2/16=3) IOP(4/21=200)";//7/9/09: similar to current number
				 *  9 with keyword small; this addresses case of VCSJVABXVCFDRA-UHFFFAOYAI 
				 *  (InChI=1/C8H19O5Si/c1-5-10-8(4)13-14(9,11-6-2)12-7-3/h8H,5-7H2,1-4H3)
				 *  
				 *  //used for troublesome C5H7J case; note that before fixing, 
				 *  I got errors like the following: "Incomplete coordinate system.  
				 *  Try restarting with Geom=Check Guess=Read Opt=(ReadFC,NewRedundant) 
				 *  Incomplete coordinate system. Error termination via Lnk1e in l103.exe"; 
				 *  we could try to restart, but it is probably preferrable to have each keyword 
				 *  combination standalone; another keyword that may be helpful if additional problematic 
				 *  cases are encountered is opt=small; 6/9/09 note: originally,
				 *   this had # pm3 opt=(verytight,gdiis,calcall) freq IOP(2/16=3)" 
				 *   (with freq keyword), but I discovered that in this case, there are 
				 *   two thermochemistry sections and cclib parses frequencies twice, 
				 *   giving twice the number of desired frequencies and hence produces 
				 *   incorrect thermo; this turned up on C5H6JJ isomer
				 */
				inpKeyStr+="# pm3 opt=(verytight,gdiis,calcall) IOP(2/16=3)";
				break;
			case 16:
				/**
				 * gmagoon 7/3/09: it is probably best to retire this keyword combination in
				 *  light of the similar combination below         //else if(attemptNumber==6)
				 *   inpKeyStr+="# pm3 opt=(verytight,gdiis,calcall,small) IOP(2/16=3) 
				 *   IOP(4/21=2)";//6/10/09: worked for OJZYSFFHCAPVGA-UHFFFAOYAK 
				 *   (InChI=1/C5H7/c1-3-5-4-2/h1,4H2,2H3) case; IOP(4/21) keyword was key
				 *   
				 *   //6/29/09: worked for troublesome ketene case: 
				 *   CCGKOQOJPYTBIH-UHFFFAOYAO (InChI=1/C2H2O/c1-2-3/h1H2) 
				 *   (could just increase number of iterations for similar keyword 
				 *   combination above (#6 at the time of this writing), allowing symmetry,
				 *    but nosymm seemed to reduce # of iterations; I think one of nosymm or 
				 *    higher number of iterations would allow the similar keyword combination 
				 *    to converge; both are included here for robustness)
				 */
				inpKeyStr+="# pm3 opt=(verytight,gdiis,calcall,small,maxcyc=200) IOP(2/16=3) IOP(4/21=2) nosymm";
				break;
			case 17:
				/**
				 * //7/1/09: added for case of ZWMVZWMBTVHPBS-UHFFFAOYAEmult3
				 *  (InChI=1/C4H4O2/c1-3-5-6-4-2/h1-2H2/mult3)
				 */
				inpKeyStr+="# pm3 opt=(verytight,gdiis,calcall,small) IOP(2/16=3) nosymm";
				break;
			case 0:
				/**
				 * //6/10/09: used to address troublesome 
				 * FILUFGAZMJGNEN-UHFFFAOYAImult3 case (InChI=1/C5H6/c1-3-5-4-2/h3H,1H2,2H3/mult3)
				 */
				inpKeyStr+="# pm3 opt=(calcall,small,maxcyc=100) IOP(2/16=3)"; 
				break;

			default:
				throw new Exception();//this point should not be reached

			}

			/**
			 * if(multiplicity == 3) inpKeyStr+= " guess=mix"; 
			 * //assumed to be triplet biradical...use guess=mix to perform unrestricted ;
			 *  nevermind...I think this would only be for singlet biradicals based on 
			 *  http://www.gaussian.com/g_tech/g_ur/k_guess.htm
			 */
			if (QMTP.usePolar) inpKeyStr += " polar";
			FileWriter fw = new FileWriter(inputKeywords);
			fw.write(inpKeyStr);
			fw.close();
		}
		catch(Exception e){
			String err = "Error in writing inputkeywords.txt \n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}
		
		keywords = new HashMap<String, String>();
		keywords.put(G03PM3KEYWORDS.INPUT, inpKeyStr);
		
		return keywords;
	}


}
