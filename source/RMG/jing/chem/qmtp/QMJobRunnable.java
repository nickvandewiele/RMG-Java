package jing.chem.qmtp;

/**
 * Interface to which all quantum chemistry packages should obey.
 * @author nmvdewie
 *
 */
public interface QMJobRunnable {

	/**
	 * run the quantum chemistry package by initiating a Process and pointing to 
	 * the external executable
	 * @return
	 */
	public int run();
	
	/**
	 * Check whether the Process has terminated succesfully
	 * @param job
	 * @return
	 */
	public int check(Process job);
}
