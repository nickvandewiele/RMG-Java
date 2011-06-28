package jing.chem.qmtp;

/**
 * supertype for all quantum chemistry packages
 * @author nmvdewie
 *
 */
public class QMJob {

	String name;
	String dir;
	
	String inputFileExtension;
	String outputFileExtension;
	
	//the keywords denoting the executable:
	String executable;
	
	//the command line command:
	String command;
	
	public static String qmfolder= "QMfiles/";
	
	public QMJob(){
		
	}
	public QMJob(String name, String directory){
		this.name = name;
		this.dir = directory;
	}
}
