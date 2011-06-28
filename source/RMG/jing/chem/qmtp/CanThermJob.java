package jing.chem.qmtp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStream;
import java.io.InputStreamReader;

import jing.rxnSys.Logger;

public class CanThermJob extends QMJob implements QMJobRunnable {

	public CanThermJob(String name, String directory) {
		super(name, directory);

		inputFileExtension = ".can";
		outputFileExtension = "";

		executable = "python ";

		command =executable 
		+ System.getenv("RMG")+"/source/CanTherm/source/CanTherm.py "
		+ name + inputFileExtension;

	}

	@Override
	public int run() {
		Process job = null;
		try{
			job = Runtime.getRuntime().exec(command, null, new File(qmfolder));	
		}
		catch(Exception e){
			String err = "Error in running CanTherm process \n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}
		return check(job);
	}

	@Override
	public int check(Process job) {
		//0 is arbitrary...
		return 0;
	}

}
