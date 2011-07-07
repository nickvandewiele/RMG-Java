package jing.chem.qmtp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStream;
import java.io.InputStreamReader;

import jing.chem.QMTP;
import jing.rxnSys.Logger;

public class MM4HRJob extends QMJob implements QMJobRunnable {

	public MM4HRJob(String name, String directory) {
		super(name, directory);
		

		inputFileExtension = ".com";
		outputFileExtension = "";

		executable = "csh ";

		command = executable;
		command = command.concat(name);

	}

	@Override
	public int run() {
		Process job = null;
		try{
			job = Runtime.getRuntime().exec(command, null, new File(qmfolder));
		}
		catch(Exception e){
			String err = "Error in running MM4HR process \n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}
		return check(job);
	}

	@Override
	public int check(Process job) {
		int successFlag=0;
		try{
			//check for errors and display the error if there is one
			//	    InputStream is = mm4Proc.getErrorStream();
			//	    InputStreamReader isr = new InputStreamReader(is);
			//	    BufferedReader br = new BufferedReader(isr);
			//	    String line=null;
			//	    while ( (line = br.readLine()) != null) {
			//		line = line.trim();
			//		if(!line.equals("STOP   statement executed")){//string listed here seems to be typical
			//		    Logger.error(line);
			//		    flag=1;
			//		}
			//	    }
			InputStream is = job.getInputStream();
			InputStreamReader isr = new InputStreamReader(is);
			BufferedReader br = new BufferedReader(isr);
			String line=null;
			while ( (line = br.readLine()) != null) {
				//do nothing
			}
			//if there was an error, indicate that an error was obtained
			//	    if(flag==1){
			//		Logger.info("MM4 process received error (see above) on " + name);
			//	    }


			int exitValue = job.waitFor();
			job.getErrorStream().close();
			job.getOutputStream().close();
			br.close();
			isr.close();
			is.close();
		}
		catch(Exception e){
			String err = "Error in running MM4 rotor process \n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}
		
		/**
		 * TODO failure in MM4HR Job should be checked, is not done up to now!!!
		 */
		return successFlag = 1;
	}

}
