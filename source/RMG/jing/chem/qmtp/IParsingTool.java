package jing.chem.qmtp;

import jing.chem.ChemGraph;

/**
 * This interface is for all tools that are able to interpret computational
 * chemistry output files like .log from Gaussian.
 * 
 * The information is stored in a {@link IQMData} object and is passed to a 
 * {@link QMParser} type that serves as a 'mailman' object. The information in the {@link IQMData} type
 * is then fed to a {@link TDPropertiesCalculator} type that calculates the variables we are interested
 * in.
 * @author nmvdewie
 *
 */
public interface IParsingTool {

	public IQMData parse(Process job, ChemGraph p_chemGraph);

}
