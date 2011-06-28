package jing.chem.qmtp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStream;
import java.io.InputStreamReader;

import jing.chem.QMTP;
import jing.rxnSys.Logger;

public class MOPACJob extends QMJob implements QMJobRunnable {

	public MOPACJob(String name, String directory) {
		super(name, directory);

		inputFileExtension = ".mop ";
		outputFileExtension = ".out";

		executable = System.getenv("MOPAC_LICENSE")+"MOPAC2009.exe ";

		command = executable;
		//specify the input file; space is important
		command=command.concat(dir+"/"+name+inputFileExtension);
		//specify the output file
		command=command.concat(dir+"/"+name+".out");

	}

	@Override
	public int run() {
		Process job = null;
		try{
			job = Runtime.getRuntime().exec(command);
		}
		catch(Exception e){
			String err = "Error in running MOPAC process \n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}
		return check(job);
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
				Logger.info("MOPAC process received error (see above) on " + name);
			}
			int exitValue = job.waitFor();
			job.getInputStream().close();
			job.getOutputStream().close();
			br.close();
			isr.close();
			is.close();
		}
		catch(Exception e){
			String err = "Error in running MOPAC process \n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}
		//look in the output file to check for the successful termination of the calculation (this is a trimmed down version of what appears in successfulMOPACResultExistsQ (it doesn't have the InChI check)
		File file = new File(dir+"/"+name+".out");
		int failureFlag=1;//flag (1 or 0) indicating whether the MOPAC job failed
		int failureOverrideFlag=0;//flag (1 or 0) to override success as measured by failureFlag
		if(file.exists()){//if the file exists, do further checks; otherwise, we will skip to final statement and return false
			try{
				FileReader in = new FileReader(file);
				BufferedReader reader = new BufferedReader(in);
				String line=reader.readLine();
				while(line!=null){
					String trimLine = line.trim();
					if (trimLine.equals("DESCRIPTION OF VIBRATIONS")){//check for this line; if it is here, check for negative frequencies
						//if(!MopacFileContainsNegativeFreqsQ(name, directory)) failureFlag=0;
						failureFlag=0;
					}
					//negative frequencies notice example:
					//         NOTE: SYSTEM IS NOT A GROUND STATE, THEREFORE ZERO POINT
					//         ENERGY IS NOT MEANINGFULL. ZERO POINT ENERGY PRINTED
					//         DOES NOT INCLUDE THE  2 IMAGINARY FREQUENCIES
					else if (trimLine.endsWith("IMAGINARY FREQUENCIES")){
						Logger.info("*****Imaginary freqencies found:");
						failureOverrideFlag=1;
					}
					else if (trimLine.equals("EXCESS NUMBER OF OPTIMIZATION CYCLES")){//exceeding max cycles error
						failureOverrideFlag=1;
					}
					else if (trimLine.equals("NOT ENOUGH TIME FOR ANOTHER CYCLE")){//timeout error
						failureOverrideFlag=1;
					}
					line=reader.readLine();
				}
				reader.close();
				in.close();
			}
			catch(Exception e){
				String err = "Error in reading MOPAC output file \n";
				err += e.toString();
				Logger.logStackTrace(e);
				System.exit(0);
			}
		}
		if(failureOverrideFlag==1) failureFlag=1; //job will be considered a failure if there are imaginary frequencies or if job terminates to to excess time/cycles
		//if the failure flag is 0 and there are no negative frequencies, the process should have been successful
		if (failureFlag==0) successFlag=1;

		return successFlag;
	}

}
