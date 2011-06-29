package jing.chem.qmtp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

import jing.chem.ChemGraph;
import jing.chem.ThermoData;
import jing.rxnSys.Logger;

/**
 * The CanTherm parser is deviating from the other parsers' architectures
 * because it directly produces 
 * thermodynamic properties without the need to actually
 * calculate the properties based on computational chemistry data
 * @author nmvdewie
 *
 */
public class CanThermParser {

	String name;
	
	String directory;
	
	ChemGraph p_chemGraph;
	
	String inputFileExtension = ".canout";
	
	public CanThermParser(String name, String directory, ChemGraph p_chemGraph){
		
		this.name = name;
		
		this.directory = directory;
		
		this.p_chemGraph = p_chemGraph;
		
	}
	
	public ThermoData parse(){
		//5. read CanTherm output
		Double Hf298 = null;
		Double S298 = null;
		Double Cp300 = null;
		Double Cp400 = null;
		Double Cp500 = null;
		Double Cp600 = null;
		Double Cp800 = null;
		Double Cp1000 = null;
		Double Cp1500 = null;
		File file = new File(directory+"/"+name+inputFileExtension);
		try{
			FileReader in = new FileReader(file);
			BufferedReader reader = new BufferedReader(in);
			String line=reader.readLine();
			while(!line.startsWith("Hf298 S298 Cps:")){//get to the end of the file with the data we want
				line=reader.readLine();
			}
			String[] split = reader.readLine().trim().split("\\s+");//read the next line, which contains the information we want
			Hf298 = Double.parseDouble(split[0]);
			S298 = Double.parseDouble(split[1]);
			Cp300 = Double.parseDouble(split[2]);
			Cp400 = Double.parseDouble(split[3]);
			Cp500 = Double.parseDouble(split[4]);
			Cp600 = Double.parseDouble(split[5]);
			Cp800 = Double.parseDouble(split[7]);
			Cp1000 = Double.parseDouble(split[9]);
			Cp1500 = Double.parseDouble(split[14]);
			reader.close();
			in.close();
		}
		catch(Exception e){
			String err = "Error in reading CanTherm .canout file \n";
			err += e.toString();
			Logger.logStackTrace(e);
			System.exit(0);
		}
		
		ThermoData result = new ThermoData(Hf298,S298,Cp300,Cp400,Cp500,Cp600,Cp800,Cp1000,Cp1500,3,1,1,"MM4 calculation; includes CanTherm analysis of force-constant matrix");//this includes rough estimates of uncertainty
		result.setSource("MM4 calculation with CanTherm analysis");
		Logger.info("Thermo for " + name + ": "+ result.toString());//print result, at least for debugging purposes

		return result;
	}
}
