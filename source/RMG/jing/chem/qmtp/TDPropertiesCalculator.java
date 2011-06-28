package jing.chem.qmtp;


import jing.chem.ThermoData;
import jing.rxnSys.Logger;

/**
 * this calculator processes information on chemical species from a computation methods
 * like quantum chemistry methods and converts it into information on thermodynamic properties 
 * of this species, which will later be used for kinetic modelling purposes.
 * @author nmvdewie
 *
 */
public class TDPropertiesCalculator {

	IQMData data;

	PointGroup pointGroup;

	public TDPropertiesCalculator(IQMData data){
		this.data = data;
	}
	public void determinePointGroup(){
		//determine point group using the SYMMETRY Program

		PointGroupCalculator pgc = new PointGroupCalculator(data);
		pointGroup = pgc.calculate();
	}
	/**
	 * gets the statistical correction for S in dimensionless units (divided by R)
	 * @return
	 */
	public Double calculateSymmetryCorrection(){
		//determine statistical correction factor for 1. external rotational symmetry (affects rotational partition function) and 2. chirality (will add R*ln2 to entropy) based on point group
		//ref: http://cccbdb.nist.gov/thermo.asp
		//assumptions below for Sn, T, Th, O, I seem to be in line with expectations based on order reported at: http://en.wikipedia.org/w/index.php?title=List_of_character_tables_for_chemically_important_3D_point_groups&oldid=287261611 (assuming order = symmetry number * 2 (/2 if chiral))...this appears to be true for all point groups I "know" to be correct
		//minor concern: does SYMMETRY appropriately calculate all Sn groups considering 2007 discovery of previous errors in character tables (cf. Wikipedia article above)
		return QMConstants.R*Math.log(pointGroup.symmetryNumber);
	}

	public Double calculateChiralityCorrection(){
		if(pointGroup.chiral){
			return QMConstants.R*Math.log(2);
		}
		else return 0.;
	}

	/**
	 * 		//calculate thermo quantities using stat. mech. equations
	 * @return
	 */
	public ThermoData calculate() {

		//we will use number of atoms from above (alternatively, we could use the chemGraph); this is needed to test whether the species is monoatomic
		double Hf298, S298, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500;

		Hf298 = data.getEnergy()*QMConstants.Hartree_to_kcal;

		/*
		 * electronic + translation; note use of 10^5 Pa for standard pressure; 
		 * also note that molecular mass needs to be divided by 1000 for kg units
		 */
		S298 = calcElecS() + calcTransS(298.0);

		Cp300 = calcTransCp();
		Cp400 = calcTransCp();
		Cp500 = calcTransCp();
		Cp600 = calcTransCp();
		Cp800 = calcTransCp();
		Cp1000 = calcTransCp();
		Cp1500 = calcTransCp();

		/*
		 * //include statistical correction and rotational
		 *  (without symmetry number, vibrational contributions if species is polyatomic
		 */
		if(data.getNumberOfAtoms()>1){

			Cp300 += calcRotCp() + calcVibCp(300.);
			Cp400 += calcRotCp() + calcVibCp(400.); 
			Cp500 += calcRotCp() + calcVibCp(500.); 
			Cp600 += calcRotCp() + calcVibCp(600.); 
			Cp800 += calcRotCp() + calcVibCp(800.); 
			Cp1000 += calcRotCp() + calcVibCp(1000.); 
			Cp1500 += calcRotCp() + calcVibCp(1500.); 

			S298 = S298 
			-calculateSymmetryCorrection()
			+calculateChiralityCorrection()
			+calcRotS(298.15)
			+calcVibS(298.15);
		}
		return new ThermoData(Hf298,S298,Cp300,Cp400,Cp500,Cp600,Cp800,Cp1000,Cp1500,5,1,1,"PM3 or MM4 calculation");//this includes rough estimates of uncertainty
	}
	/**
	 * electronic partition function
	 * @return
	 */
	public double calcElecS() {
		return QMConstants.R * Math.log(data.getGroundStateDegeneracy());
	}
	//gmagoon 6/8/09
	//calculate the vibrational contribution (divided by R, dimensionless) at temperature, T, in Kelvin to entropy
	//p_freqs in cm^-1; c in cm/s; k in J/K; h in J-s
	//ref.: http://cccbdb.nist.gov/thermo.asp
	public double calcVibS(double temperature){
		double Scontrib = 0;
		double dr;
		for(int i=0; i < data.getFrequencies().size(); i++){
			double freq = (Double)data.getFrequencies().get(i);
			dr = QMConstants.h*QMConstants.c*freq/(QMConstants.k*temperature); //frequently used dimensionless ratio
			Scontrib = Scontrib - Math.log(1.-Math.exp(-dr))+dr*Math.exp(-dr)/(1.-Math.exp(-dr));
		}

		return QMConstants.R*Scontrib;
	}

	//gmagoon 6/8/09
	//calculate the vibrational contribution (divided by R, dimensionless) at temperature, T, in Kelvin to heat capacity, Cp
	//p_freqs in cm^-1; c in cm/s; k in J/K; h in J-s
	//ref.: http://cccbdb.nist.gov/thermo.asp
	public double calcVibCp(double temperature){
		double Cpcontrib = 0;
		double dr;
		for(int i=0; i < data.getFrequencies().size(); i++){
			double freq = (Double)data.getFrequencies().get(i);
			dr = QMConstants.h*QMConstants.c*freq/(QMConstants.k*temperature); //frequently used dimensionless ratio
			Cpcontrib = Cpcontrib + Math.pow(dr, 2.)*Math.exp(-dr)/Math.pow(1.-Math.exp(-dr),2.);
		}

		return QMConstants.R*Cpcontrib;
	}

	//gmagoon 6/23/10
	//calculate the vibrational contribution (divided by R, units of K) at temperature, T, in Kelvin to Hthermal (ZPE not included)
	//p_freqs in cm^-1; c in cm/s; k in J/K; h in J-s
	//we need to ignore zero frequencies, as MM4 does, to avoid dividing by zero; based on equation for Ev in http://www.gaussian.com/g_whitepap/thermo.htm; however, we ignore zero point contribution to be consistent with the input that CanTherm currently takes; note, however, that this could introduce some small inaccuracies, as the frequencies may be slightly different in CanTherm vs. MM4, particularly for a frequency < 7.7 cm^-1 (treated as zero in MM4)
	public double calcVibH(double temperature){
		double Hcontrib = 0;
		double dr;
		for(int i=0; i < data.getFrequencies().size(); i++){
			double freq = (Double)data.getFrequencies().get(i);
			if(freq > 0.0){//ignore zero frequencies as MM4 does
				dr = QMConstants.h*QMConstants.c*freq/(QMConstants.k*temperature); //frequently used dimensionless ratio
				Hcontrib = Hcontrib + dr*temperature/(Math.exp(dr) - 1.);
			}
		}

		return Hcontrib;
	}
	/**
	 * Translational contribution to the molar entropy at standard conditions (1 bar)
	 * 
	 * @param temperature
	 * @return
	 */
	public double calcTransS(double temperature){
		return QMConstants.R * (3./2.*Math.log(2.*Math.PI * data.getMolecularMass()
				/(1000.*QMConstants.Na*Math.pow(QMConstants.h,2.)))
				+5./2.*Math.log(QMConstants.k*temperature)
				-Math.log(100000.)+5./2.);
	}

	/**
	 * translational contribution to standard heat capacity at constant pressure: 
	 * @return
	 */
	public double calcTransCp(){
		return 5./2.*QMConstants.R;
	}

	public double calcRotCp(){
		if(pointGroup.isLinear()){
			return QMConstants.R;
		}
		else return 3./2. * QMConstants.R;
	}

	/**
	 * rotational contribution to entropy
	 * 
	 * for linear molecules, one rotational constant will be zero, 
	 * the other two will be identical
	 * @return
	 */
	public double calcRotS(double temperature){
		double rotCons;
		if(pointGroup.isLinear()){
			if(data.getRotationalConstants()[0] > 0.0001) rotCons = data.getRotationalConstants()[0];
			else rotCons = data.getRotationalConstants()[1];

			return QMConstants.R*(Math.log(QMConstants.k*temperature/(QMConstants.h*rotCons))+1);
		}
		else{
			return QMConstants.R*(3./2.*Math.log(QMConstants.k*temperature/QMConstants.h)
					-1./2.*Math.log(data.getRotationalConstants()[0]*data.getRotationalConstants()[1]*data.getRotationalConstants()[2]/Math.PI)+3./2.)
					+QMConstants.R*calcVibS(temperature);
		}

	}
}
