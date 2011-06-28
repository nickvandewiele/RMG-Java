package jing.chem.qmtp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStream;
import java.io.InputStreamReader;

import jing.chem.QMTP;
import jing.rxnSys.Logger;

public class MM4Job extends QMJob implements QMJobRunnable {

	public MM4Job(String name, String directory) {
		super(name, directory);

		inputFileExtension = ".com";
		outputFileExtension = ".mm4out";

		executable = "csh ";
		// File script = new File(qmfolder+command);
		// command = "./"+command;
		// script.setExecutable(true);
		command = executable;
		command = command.concat(name+inputFileExtension);

	}

	@Override
	public int run() {
		Process job = null;
		try{
			job = Runtime.getRuntime().exec(command, null, new File(qmfolder));
		}
		catch(Exception e){
			String err = "Error in running Gaussian process \n";
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
			String err = "Error in running MM4 process \n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}
		//look in the output file to check for the successful termination of the MM4 calculation (cf. successfulMM4ResultExistsQ)
		File file = new File(dir+"/"+name+outputFileExtension);
		int failureFlag=1;//flag (1 or 0) indicating whether the MM4 job failed
		int failureOverrideFlag=0;//flag (1 or 0) to override success as measured by failureFlag
		if(file.exists()){//if the file exists, do further checks; otherwise, we will skip to final statement and return false
			try{
				FileReader in = new FileReader(file);
				BufferedReader reader = new BufferedReader(in);
				String line=reader.readLine();
				while(line!=null){
					String trimLine= line.trim();
					if (trimLine.equals("STATISTICAL THERMODYNAMICS ANALYSIS")){
						failureFlag = 0;
					}
					else if (trimLine.endsWith("imaginary frequencies,")){//read the number of imaginary frequencies and make sure it is zero
						String[] split = trimLine.split("\\s+");
						if (Integer.parseInt(split[3])>0){
							Logger.info("*****Imaginary freqencies found:");
							failureOverrideFlag=1;
						}
					}
					else if (trimLine.contains("             0.0     (fir )")){
						if (QMTP.useCanTherm){//zero frequencies are only acceptable when CanTherm is used
							Logger.info("*****Warning: zero freqencies found (values lower than 7.7 cm^-1 are rounded to zero in MM4 output); CanTherm should hopefully correct this:");
						}
						else{
							Logger.info("*****Zero freqencies found:");
							failureOverrideFlag=1;
						}
					}
					line=reader.readLine();
				}
				reader.close();
				in.close();
			}
			catch(Exception e){
				String err = "Error in reading MM4 output file \n";
				err += e.toString();
				Logger.logStackTrace(e);
				System.exit(0);
			}
		}
		//if the failure flag is still 0, the process should have been successful
		if(failureOverrideFlag==1) failureFlag=1; //job will be considered a failure if there are imaginary frequencies or if job terminates to to excess time/cycles
		//if the failure flag is 0 and there are no negative frequencies, the process should have been successful
		if (failureFlag==0) successFlag=1;

		return successFlag;
	}

}
