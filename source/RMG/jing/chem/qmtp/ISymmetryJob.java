package jing.chem.qmtp;

/**
 * Interface to all algorithms that perceive symmetry numbers based on 
 * 3D coords of a molecule
 * @author nmvdewie
 *
 */
public interface ISymmetryJob {

	public PointGroup calculate();
}
