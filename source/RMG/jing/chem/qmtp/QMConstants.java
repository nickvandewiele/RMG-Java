package jing.chem.qmtp;

/**
 * Container for static final quantum chemistry constants
 * @author Nick Vandewiele
 *
 */
public class QMConstants {

	/**
	 * ideal gas constant in cal/mol-K 
	 * (does this appear elsewhere in RMG, so I don't need to reuse it?)
	 */
	public static final double R = 1.9872;
	/**
	 * Avagadro's number; 
	 * cf. http://physics.nist.gov/cgi-bin/cuu/Value?na|search_for=physchem_in!
	 */
	public static final double Na = 6.02214179E23; 
	public static double Hartree_to_kcal = 627.5095; //conversion from Hartree to kcal/mol taken from Gaussian thermo white paper
	public static double k = 1.3806504E-23;//Boltzmann's constant in J/K; cf. http://physics.nist.gov/cgi-bin/cuu/Value?na|search_for=physchem_in!
	public static double h = 6.62606896E-34;//Planck's constant in J-s; cf. http://physics.nist.gov/cgi-bin/cuu/Value?h|search_for=universal_in!
	public static double c = 299792458. *100;//speed of light in vacuum in cm/s, cf. http://physics.nist.gov/cgi-bin/cuu/Value?c|search_for=universal_in!
	public final static double ENTHALPY_HYDROGEN = 52.1; //needed for HBI

}
