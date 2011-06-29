package jing.chem.qmtp;

import jing.chem.ChemGraph;
import jing.chem.ThermoData;
import jing.rxnSys.Logger;

/**
 * Supertype for all quantum chemistry parsers that retrieve
 * return a {@link ThermoData} object with TD information from a computational
 * chemistry output file.
 * 
 * The actual parsing is 'outsourced' to parser tools, like CCLib, CanTherm.
 * 
 * Therefore, an attribute {@link ParsingTool} contains an implementation like
 * CCLib, CanTherm, that will be specified in the constructor.
 * @author nmvdewie
 *
 */
public abstract class QMParser implements QMParsable{

	String name;
	
	String dir;
	
	ChemGraph p_chemGraph;
	
	String inputFileExtension;
	
	//the keywords denoting the executable:
	String executable;
	
	//the command line command:
	String command;
	
	//directory of parsing scripts:
	String scripts = "/scripts/";
	
	//parsing script file:
	String scriptFile;
	
	IParsingTool parsingTool;
	
	IQMData data;
	
	public QMParser(String name, String directory, ChemGraph p_chemGraph){
		this.name = name;
		
		this.dir = directory;
		
		this.p_chemGraph = p_chemGraph;
	}
	@Override
	public ThermoData parse() {
	
		//        //parse the Gaussian file using cclib
		//        int natoms = 0; //number of atoms from Gaussian file; in principle, this should agree with number of chemGraph atoms
		//        ArrayList atomicNumber = new ArrayList(); //vector of atomic numbers (integers) (apparently Vector is thread-safe; cf. http://answers.yahoo.com/question/index?qid=20081214065127AArZDT3; ...should I be using this instead?)
		//        ArrayList x_coor = new ArrayList(); //vectors of x-, y-, and z-coordinates (doubles) (Angstroms) (in order corresponding to above atomic numbers)
		//        ArrayList y_coor = new ArrayList();
		//        ArrayList z_coor = new ArrayList();
		//        double energy = 0; //PM3 energy (Hf298) in Hartree
		//        double molmass = 0; //molecular mass in amu
		//        ArrayList freqs = new ArrayList(); //list of frequencies in units of cm^-1
		//        double rotCons_1 = 0;//rotational constants in (1/s)
		//        double rotCons_2 = 0;
		//        double rotCons_3 = 0; 
		//        int gdStateDegen = p_chemGraph.getRadicalNumber()+1;//calculate ground state degeneracy from the number of radicals; this should give the same result as spin multiplicity in Gaussian input file (and output file), but we do not explicitly check this (we could use "mult" which cclib reads in if we wanted to do so); also, note that this is not always correct, as there can apparently be additional spatial degeneracy for non-symmetric linear molecules like OH radical (cf. http://cccbdb.nist.gov/thermo.asp)
		//        try{   
		//            File runningdir=new File(directory);
		//            String command = "c:/Python25/python.exe c:/Python25/GaussianPM3ParsingScript.py ";//this should eventually be modified for added generality
		//            String logfilepath=directory+"/"+name+".log";
		//            command=command.concat(logfilepath);
		//            Process cclibProc = Runtime.getRuntime().exec(command, null, runningdir);
		//            //read the stdout of the process, which should contain the desired information in a particular format
		//            InputStream is = cclibProc.getInputStream();
		//            InputStreamReader isr = new InputStreamReader(is);
		//            BufferedReader br = new BufferedReader(isr);
		//            String line=null;
		//            //example output:
		////            C:\Python25>python.exe GaussianPM3ParsingScript.py TEOS.out
		////            33
		////            [ 6  6  8 14  8  6  6  8  6  6  8  6  6  1  1  1  1  1  1  1  1  1  1  1  1
		////              1  1  1  1  1  1  1  1]
		////            [[ 2.049061 -0.210375  3.133106]
		////             [ 1.654646  0.321749  1.762752]
		////             [ 0.359284 -0.110429  1.471465]
		////             [-0.201871 -0.013365 -0.12819 ]
		////             [ 0.086307  1.504918 -0.82893 ]
		////             [-0.559186  2.619928 -0.284003]
		////             [-0.180246  3.839463 -1.113029]
		////             [ 0.523347 -1.188305 -1.112765]
		////             [ 1.857584 -1.018167 -1.495088]
		////             [ 2.375559 -2.344392 -2.033403]
		////             [-1.870397 -0.297297 -0.075427]
		////             [-2.313824 -1.571765  0.300245]
		////             [-3.83427  -1.535927  0.372171]
		////             [ 1.360346  0.128852  3.917699]
		////             [ 2.053945 -1.307678  3.160474]
		////             [ 3.055397  0.133647  3.403037]
		////             [ 1.677262  1.430072  1.750899]
		////             [ 2.372265 -0.029237  0.985204]
		////             [-0.245956  2.754188  0.771433]
		////             [-1.656897  2.472855 -0.287156]
		////             [-0.664186  4.739148 -0.712606]
		////             [-0.489413  3.734366 -2.161038]
		////             [ 0.903055  4.016867 -1.112198]
		////             [ 1.919521 -0.229395 -2.269681]
		////             [ 2.474031 -0.680069 -0.629949]
		////             [ 2.344478 -3.136247 -1.273862]
		////             [ 1.786854 -2.695974 -2.890647]
		////             [ 3.41648  -2.242409 -2.365094]
		////             [-1.884889 -1.858617  1.28054 ]
		////             [-1.976206 -2.322432 -0.440995]
		////             [-4.284706 -1.26469  -0.591463]
		////             [-4.225999 -2.520759  0.656131]
		////             [-4.193468 -0.809557  1.112677]]
		////            -14.1664924726
		////            [    9.9615    18.102     27.0569    31.8459    39.0096    55.0091
		////                66.4992    80.4552    86.4912   123.3551   141.6058   155.5448
		////               159.4747   167.0013   178.5676   207.3738   237.3201   255.3487
		////               264.5649   292.867    309.4248   344.6503   434.8231   470.2074
		////               488.9717   749.1722   834.257    834.6594   837.7292   839.6352
		////               887.9767   892.9538   899.5374   992.1851  1020.6164  1020.8671
		////              1028.3897  1046.7945  1049.1768  1059.4704  1065.1505  1107.4001
		////              1108.1567  1109.0466  1112.6677  1122.7785  1124.4315  1128.4163
		////              1153.3438  1167.6705  1170.9627  1174.9613  1232.1826  1331.8459
		////              1335.3932  1335.8677  1343.9556  1371.37    1372.8127  1375.5428
		////              1396.0344  1402.4082  1402.7554  1403.2463  1403.396   1411.6946
		////              1412.2456  1412.3519  1414.5982  1415.3613  1415.5698  1415.7993
		////              1418.5409  2870.7446  2905.3132  2907.0361  2914.1662  2949.2646
		////              2965.825   2967.7667  2971.5223  3086.3849  3086.3878  3086.6448
		////              3086.687   3089.2274  3089.4105  3089.4743  3089.5841  3186.0753
		////              3186.1375  3186.3511  3186.365 ]
		////            [ 0.52729  0.49992  0.42466]
		////note: above example has since been updated to print molecular mass; also frequency and atomic number format has been updated
		//            String [] stringArray;
		//            natoms = Integer.parseInt(br.readLine());//read line 1: number of atoms
		//            stringArray = br.readLine().replace("[", "").replace("]","").trim().split(",\\s+");//read line 2: the atomic numbers (first removing braces)
		//           // line = br.readLine().replace("[", "").replace("]","");//read line 2: the atomic numbers (first removing braces)
		//           // StringTokenizer st = new StringTokenizer(line); //apprently the stringTokenizer class is deprecated, but I am having trouble getting the regular expressions to work properly
		//            for(int i=0; i < natoms; i++){
		//               // atomicNumber.add(i,Integer.parseInt(stringArray[i]));
		//                atomicNumber.add(i,Integer.parseInt(stringArray[i]));
		//            }
		//            for(int i=0; i < natoms; i++){
		//                stringArray = br.readLine().replace("[", "").replace("]","").trim().split("\\s+");//read line 3+: coordinates for atom i; used /s+ for split; using spaces with default limit of 0 was giving empty string
		//                x_coor.add(i,Double.parseDouble(stringArray[0]));
		//                y_coor.add(i,Double.parseDouble(stringArray[1]));
		//                z_coor.add(i,Double.parseDouble(stringArray[2]));
		//            }
		//            energy = Double.parseDouble(br.readLine());//read next line: energy
		//            molmass = Double.parseDouble(br.readLine());//read next line: molecular mass (in amu)
		//            if (natoms>1){//read additional info for non-monoatomic species
		//                stringArray = br.readLine().replace("[", "").replace("]","").trim().split(",\\s+");//read next line: frequencies
		//                for(int i=0; i < stringArray.length; i++){
		//                    freqs.add(i,Double.parseDouble(stringArray[i]));
		//                }
		//                stringArray = br.readLine().replace("[", "").replace("]","").trim().split("\\s+");//read next line rotational constants (converting from GHz to Hz in the process)
		//                rotCons_1 = Double.parseDouble(stringArray[0])*1000000000;
		//                rotCons_2 = Double.parseDouble(stringArray[1])*1000000000;
		//                rotCons_3 = Double.parseDouble(stringArray[2])*1000000000;
		//            }
		//            while ( (line = br.readLine()) != null) {
		//                //do nothing (there shouldn't be any more information, but this is included to get all the output)
		//            }
		//            int exitValue = cclibProc.waitFor();
		//        }
		//        catch (Exception e) {
		//				Logger.logStackTrace(e);
		//            String err = "Error in running ccLib Python process \n";
		//            err += e.toString();
		//            Logger.critical(err);
		//            System.exit(0);
		//        } 
		//   
		//        ThermoData result = calculateThermoFromPM3Calc(natoms, atomicNumber, x_coor, y_coor, z_coor, energy, molmass, freqs, rotCons_1, rotCons_2, rotCons_3, gdStateDegen);
		//        Logger.info("Thermo for " + name + ": "+ result.toString());//print result, at least for debugging purposes
		//        return result;

		Process parserJob = run();
		
		read(parserJob);
		
		TDPropertiesCalculator calculator = new TDPropertiesCalculator(data);
		
		return calculator.calculate();
		
	}

	@Override
	public Process run() {
		Process job = null;
		try{
			job = Runtime.getRuntime().exec(command);	
		}

		catch (Exception e) {
			Logger.logStackTrace(e);
			String err = "Error in running ccLib Python process \n";
			err += e.toString();
			Logger.critical(err);
			System.exit(0);
		} 
		return job;

	}

	@Override
	public IQMData read(Process job) {
		data = parsingTool.parse(job, p_chemGraph);	
		return data;

	}
	public IQMData getQMData() {
		return data;
	}
}
