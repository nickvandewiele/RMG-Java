package jing.chem.qmtp;

/**
 * Wrapper type to determine molecular symmetry point groups based on 3D coords information.
 * 
 * Will point to a specific algorithm, like SYMMETRY that is able to do this.
 * @author nmvdewie
 *
 */
public class PointGroupCalculator {

	ISymmetryJob calculator;
	
	IQMData data;
	
	public  PointGroup calculate() {
		return calculator.calculate();
	}
	
	public PointGroupCalculator(IQMData qmdata){
		
		calculator = new SymmetryJob(qmdata);
		
		this.data = qmdata;
	}

}
