package jing.chem.qmtp;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.InputStream;
import java.io.InputStreamReader;

import jing.rxnSys.Logger;

public class GaussianJob extends QMJob implements QMJobRunnable {

	public GaussianJob(String name, String directory) {
		super(name, directory);

		inputFileExtension = ".gjf";
		outputFileExtension = ".log";

		executable = "g03 ";

		command = executable;
		command = command.concat(qmfolder+"/"+name+".gjf ");//specify the input file; space is important
		command=command.concat(qmfolder+"/"+name+".log");//specify the output file


	}

	@Override
	public int run() {
		Process gaussianProc = null;
		try{
			gaussianProc = Runtime.getRuntime().exec(command);	
		}
		catch(Exception e){
			String err = "Error in running Gaussian process \n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}
		return check(gaussianProc);
	}

	@Override
	public int check(Process job) {
		int flag = 0;
		int successFlag=0;
		try{
			//check for errors and display the error if there is one
			InputStream is = job.getErrorStream();
			InputStreamReader isr = new InputStreamReader(is);
			BufferedReader br = new BufferedReader(isr);
			String line=null;
			while ( (line = br.readLine()) != null) {
				line = line.trim();
				Logger.error(line);
				flag=1;
			}
			//if there was an error, indicate that an error was obtained
			if(flag==1){
				Logger.info("Gaussian process received error (see above) on " + name);
			}
			int exitValue = job.waitFor();
			job.getInputStream().close();
			job.getOutputStream().close();
			br.close();
			isr.close();
			is.close();
		}
		catch(Exception e){
			String err = "Error in running Gaussian process \n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}
		//look in the output file to check for the successful termination of the Gaussian calculation
		//failed jobs will contain the a line beginning with " Error termination" near the end of the file
		int failureFlag=0;
		int completeFlag = 0;
		String errorLine = "";//string to store the error
		try{
			FileReader in = new FileReader(dir+"/"+name+outputFileExtension);
			BufferedReader reader = new BufferedReader(in);
			String line=reader.readLine();
			while(line!=null){
				if (line.startsWith(" Error termination ")){
					failureFlag=1;
					errorLine = line.trim();
					Logger.info("*****Error in Gaussian log file: "+errorLine);//print the error (note that in general, I think two lines will be printed)
				}
				else if (line.startsWith(" ******")){//also look for imaginary frequencies
					if (line.contains("imaginary frequencies")){
						Logger.info("*****Imaginary freqencies found:");
						failureFlag=1;
					}
				}
				else if (line.startsWith(" Normal termination of Gaussian 03 at")){//check for completion notice at end of file
					if (reader.readLine()==null){//the above line will occur once in the middle of most files (except for monatomic species), but the following line is non-critical (e.g. "Link1:  Proceeding to internal job step number  2.") so reading over this line should not cause a problem
						completeFlag=1;
					}
				}
				line=reader.readLine();
			}
			reader.close();
			in.close();
		}
		catch(Exception e){
			String err = "Error in reading Gaussian log file \n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}
		//if the complete flag is still 0, the process did not complete and is a failure
		if (completeFlag==0) failureFlag=1;
		//if the failure flag is still 0, the process should have been successful
		if (failureFlag==0) successFlag=1;

		return successFlag;
	}

}
