package jing.chem.qmtp;

import jing.chem.ThermoData;

/**
 * Interface for readers that read quantum chemistry output files,
 * and parses the information of this output file into an RMG object
 * @author nmvdewie
 *
 */
public interface QMParsable {

	/**
	 * wrapper method that first runs the parser, and then collects the information
	 * we are interested in.
	 */
	public ThermoData parse();
	
	/**
	 * points to the external executable to interpret the 
	 * quantum chemistry output file
	 * 
	 * @return the Java Process wrapper 
	 */
	public Process run();
	/**
	 * Uses the Process and collects (reads) the information we are interested in
	 * from the parser
	 * w 
	 * @param p
	 */
	public IQMData read(Process p);
}
