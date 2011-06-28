package jing.chem.qmtp;

import java.io.File;
import java.util.Map;

/**
 * all writers of input files for quantum chemistry packages should implement
 * a writer method that writes this file.
 * @author nmvdewie
 *
 */
public interface QMInputWritable {

	/**
	 * writes the input file
	 * @return
	 */
	public File write();
	
	/**
	 * creates a map of the particular method and its associated keywords.
	 * A particular quantum chemistry method can only take one String of keywords;
	 * others may require multiple Strings of independent keywords.
	 * 
	 * The Map interface is therefore a useful container.
	 * 
	 * @return
	 */
	public Map<String, String> createKeywords();
	
	/**
	 * the method requirement to specify how the quantum chemistry input file
	 * needs to be created.
	 * @return
	 */
	public File createInputFile();
}
