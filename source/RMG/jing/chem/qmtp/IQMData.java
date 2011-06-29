package jing.chem.qmtp;

import java.util.List;
import java.util.Vector;

import quicktime.qd3d.math.Vector3D;

/**
 * This interface models a data container type that collects 
 * information on chemical properties of chemical species from
 * computation chemistry methods parser, like CCLib that is able
 * to interpret for example quantum chemistry output files of a number
 * of QM packages, like Gaussian, MOPAC.
 * 
 * This data object can then be used to calculate thermodynamic properties like
 * enthalpy, entropy, heat capacity for kinetic modelling purposes.
 * 
 * The interface demands to provide getter methods to a number of species properties
 * that are required to calculated the variables (H°, ...) we are interested in.
 * 
 * This data can then 
 * @author nmvdewie
 *
 */
public interface IQMData {
	
	public Integer getNumberOfAtoms() ;

	public List<Integer> getAtomicNumbers() ;
	
	public List<Vector3D> getThreeDCoords() ;

	public Double getEnergy() ;
	
	public Double getMolecularMass() ;

	public List<Double> getFrequencies();

	public Integer getGroundStateDegeneracy() ;

	public Double[] getRotationalConstants() ;
	
	public void setNumberOfAtoms(Integer numberOfAtoms);

	public void setAtomicNumbers(List<Integer> atomicNumbers) ;

	public void setThreeDCoords(List<Vector3D> threeDCoords) ;

	public void setEnergy(Double energy) ;

	public void setMolecularMass(Double molecularMass) ;

	public void setFrequencies(List<Double> frequencies) ;

	public void setGroundStateDegeneracy(Integer groundStateDegeneracy) ;

	public void setRotationalConstants(Double[] rotationalConstants) ;
	
	public Double getStericEnergy() ;
	public void setStericEnergy(Double stericEnergy);
}
