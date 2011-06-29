package jing.chem.qmtp;

import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;

import jing.chem.ChemGraph;
import jing.chem.QMTP;
import jing.rxnSys.Logger;

public class CanThermInputWriter extends QMInputWriter implements QMInputWritable {
	
	IQMData data;
	
	ChemGraph chemGraph;
	
	double[] dihedralMinima;
	
	Boolean forceRRHO;
	
	private class CANTHERMKEYS {
		public static final String CANTHERM = "Cantherm";
		public static final String ROTOR = "Rotor";
	};
	public CanThermInputWriter(String name, String dir){
		super(name,dir);
	}
	public CanThermInputWriter(String name, String dir, ChemGraph chemGraph, IQMData data,
			double[] dihedralMinima, boolean forceRRHO){
		this(name,dir);
		this.chemGraph = chemGraph;
		this.data = data;
		this.dihedralMinima = dihedralMinima;
		this.forceRRHO = forceRRHO;
	}

	@Override
	public File write() {
		createKeywords();

		File inputFile = createInputFile();
		
		return inputFile;
	}

	@Override
	public Map<String, String> createKeywords() {
	
		double energy = data.getEnergy();
		double stericEnergy = data.getStericEnergy();
		List freqs = data.getFrequencies();
		
		//2. compute H0/E0;  note that we will compute H0 for CanTherm by H0=Hf298(harmonicMM4)-(H298-H0)harmonicMM4, where harmonicMM4 values come from cclib parsing;  298.16 K is the standard temperature used by MM4; also note that Hthermal should include the R*T contribution (R*298.16 (enthalpy vs. energy difference)) and H0=E0 (T=0)	
		double T_MM4 = 298.16;
		//convert from Hartree to kcal/mol
		energy *= QMConstants.Hartree_to_kcal;
		//convert from Hartree to kcal/mol
		stericEnergy *= QMConstants.Hartree_to_kcal;
		//enthalpy vs. energy difference(RT)+translation(3/2*R*T) contribution to thermal energy
		double Hthermal = 5./2.*QMConstants.R*T_MM4/1000.;
		//rotational contribution
		if(chemGraph.getAtomNumber()==1) Hthermal += 0.0;
		else if (chemGraph.isLinear()) Hthermal += QMConstants.R * T_MM4/1000.;
		else Hthermal += 3./2.*QMConstants.R*T_MM4/1000;
		//vibrational contribution
		if(chemGraph.getAtomNumber()!=1)Hthermal+=QMConstants.R*calcVibH(freqs, T_MM4)/1000.;
		energy = energy - Hthermal;
		//3. write CanTherm input file
		//determine point group using the SYMMETRY Program

		PointGroupCalculator pgc = new PointGroupCalculator(data);
		PointGroup pointGroup = pgc.calculate();
		double chiralityCorrection = 0;
		if(pointGroup.chiral){
			chiralityCorrection = 2;
		}
		else {
			chiralityCorrection = 1;
		}
		double sigmaCorr = pointGroup.symmetryNumber - chiralityCorrection;
		
		String canInp = "Calculation: Thermo\n";
		canInp += "Trange: 300 100 13\n";//temperatures from 300 to 1500 in increments of 100
		canInp += "Scale: 1.0\n";//scale factor of 1
		canInp += "Mol 1\n";
		if(chemGraph.getAtomNumber()==1) canInp += "ATOM\n";
		else if (chemGraph.isLinear()) canInp+="LINEAR\n";
		else canInp+="NONLINEAR\n";
		canInp += "GEOM MM4File " + name+".mm4out\n";//geometry file; ***special MM4 treatment in CanTherm; another option would be to use mm4opt file, but CanTherm code would need to be modified accordingly
		canInp += "FORCEC MM4File "+name+".fmat\n";//force constant file; ***special MM4 treatment in CanTherm
		if (forceRRHO) name = name + "_RRHO"; //"_RRHO" will be appended to InChIKey for RRHO calcs (though we want to use unmodified name in getQMDataWithCClib above, and in GEOM and FORCEC sections
		canInp += "ENERGY "+ energy +" MM4\n";//***special MM4 treatment in CanTherm
		canInp+="EXTSYM "+sigmaCorr+"\n";//***modified treatment in CanTherm; traditional EXTSYM integer replaced by EXTSYM double, to allow fractional values that take chirality into account
		canInp+="NELEC 1\n";//multiplicity = 1; all cases we consider will be non-radicals
		int rotors = chemGraph.getInternalRotor();
		String rotInput = null;
		if(!QMTP.useHindRot || rotors==0 || forceRRHO) canInp += "ROTORS 0\n";//do not consider hindered rotors
		else{
			int rotorCount = 0;
			canInp+="ROTORS "+rotors+ " "+name+".rotinfo\n";
			canInp+="POTENTIAL separable mm4files_inertia\n";//***special MM4 treatment in canTherm;
			String rotNumbersLine=""+stericEnergy;//the line that will contain V0 (kcal/mol), and all the dihedral minima (degrees)
			rotInput = "L1: 1 2 3\n";
			LinkedHashMap rotorInfo = chemGraph.getInternalRotorInformation();
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
		
		keywords = new HashMap<String, String>();
		keywords.put(CANTHERMKEYS.CANTHERM, canInp);
		keywords.put(CANTHERMKEYS.ROTOR, rotInput);
		
		return keywords;
		
	}

	@Override
	public File createInputFile() {
		File canFile = null;
		try{
			canFile=new File(dir+"/"+name+".can");
			FileWriter fw = new FileWriter(canFile);
			fw.write(keywords.get(CANTHERMKEYS.CANTHERM));
			fw.close();
			if(keywords.get(CANTHERMKEYS.ROTOR) !=null){//write the rotor information
				File rotFile=new File(dir+"/"+name+".rotinfo");
				FileWriter fwr = new FileWriter(rotFile);
				fwr.write(CANTHERMKEYS.ROTOR);
				fwr.close();
			}
		}
		catch(Exception e){
			String err = "Error in writing CanTherm input \n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}
		
		return canFile;
	}
	//gmagoon 6/23/10
	//calculate the vibrational contribution (divided by R, units of K) at temperature, T, in Kelvin to Hthermal (ZPE not included)
	//p_freqs in cm^-1; c in cm/s; k in J/K; h in J-s
	//we need to ignore zero frequencies, as MM4 does, to avoid dividing by zero; based on equation for Ev in http://www.gaussian.com/g_whitepap/thermo.htm; however, we ignore zero point contribution to be consistent with the input that CanTherm currently takes; note, however, that this could introduce some small inaccuracies, as the frequencies may be slightly different in CanTherm vs. MM4, particularly for a frequency < 7.7 cm^-1 (treated as zero in MM4)
	public double calcVibH(List p_freqs, double p_T){
		double Hcontrib = 0;
		double dr;
		for(int i=0; i < p_freqs.size(); i++){
			double freq = (Double)p_freqs.get(i);
			if(freq > 0.0){//ignore zero frequencies as MM4 does
				dr = QMConstants.h*QMConstants.c*freq/(QMConstants.k*p_T); //frequently used dimensionless ratio
				Hcontrib = Hcontrib + dr*p_T/(Math.exp(dr) - 1.);
			}
		}

		return Hcontrib;
	}
}
