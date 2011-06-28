package jing.chem.qmtp;


/**
 * Generic RMG Exception subclassed from the top level type
 * @author nmvdewie
 *
 */
public class RMGException extends Exception {

	public RMGException(String message) {
		super(message);
	}
	
	public RMGException(String message, Throwable cause) {
        super(message, cause);
    }
}
