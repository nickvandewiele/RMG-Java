package jing.chem.qmtp;

import java.util.Map;

import jing.chem.molFile;

/**
 * supertype for all input writers for quantum chemistry methods
 * @author nmvdewie
 *
 */
public class QMInputWriter {

	public String name;
	
	public String dir;
	
	public molFile molfile;
	
	public int attemptNumber;
	
	public int multiplicity;
	
	public Map<String, String> keywords;
	
	/**
	 * for diagnostics reasons, keywords can hereby be retrieved for a specific
	 * QM Job.
	 * TODO should probably be in an interface, instead of supertype...
	 * @return
	 */
	public Map<String, String> getKeywords() {
		if (keywords != null) return keywords;
		else return null;
	}
	public QMInputWriter(String name, String directory){
		this.name = name;
		this.dir = directory;
	}
	public QMInputWriter(String name, String directory, molFile p_molfile,
			int attemptNumber, int multiplicity) {
		this(name, directory);
		this.molfile = p_molfile;
		this.attemptNumber = attemptNumber;
		this.multiplicity = multiplicity;
	}
	/**
	 * the number of keyword permutations available; update as additional options are added
	 */
	public static int scriptAttempts;
	/**
	 * we will try a second time with crude coordinates if the UFF refined coordinates do not work
	 */
	public static int maxAttemptNumber;

}
