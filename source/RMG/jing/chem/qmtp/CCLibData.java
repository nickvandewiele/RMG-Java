package jing.chem.qmtp;

import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import quicktime.qd3d.math.Vector3D;

/**
 * This data objects collects information that CCLib was able to
 * retrieve from a quantum chemistry output file
 * @author nmvdewie
 *
 */
public class CCLibData implements IQMData{
	
	/**
	 * number of atoms from Mopac file; in principle,
	 * this should agree with number of chemGraph atoms
	 */
	Integer numberOfAtoms;
	
	/**
	 * order lists of atom types and their associated atomic numbers:
	 * atomic numbers are considered as Integers 
	 */
	List<Integer> atomicNumbers;
	
	/**
	 * an ordered list of the three dimensional coordinates of the all the atoms in 
	 * the molecule.
	 * 
	 * Coordinates are floating numbers
	 */
	List<Vector3D> threeDCoords;
	
	/**
	 * Energy in Hartree
	 */
	Double energy;
	
	/**
	 * Molecular Mass in amu
	 */
	Double molecularMass; 
	
	
	/**
	 * ordered list of frequences in cm-1
	 */
	List<Double> frequencies;
	
	/**
	 * Electronic ground state degeneracy
	 * 
	 * in RMG taken as Number of radicals +1
	 */
	Integer groundStateDegeneracy;
	
	
	/**
	 * 3 Rotational constants in a fixed Double array:
	 */
	Double[] rotationalConstants = new Double[3];
	
	Double stericEnergy;
	
	public CCLibData(){
		atomicNumbers = new ArrayList<Integer>();
		
		threeDCoords = new ArrayList<Vector3D>();
		
		frequencies = new ArrayList<Double>();
		
	}

	public Integer getNumberOfAtoms() {
		return numberOfAtoms;
	}

	public void setNumberOfAtoms(Integer numberOfAtoms) {
		this.numberOfAtoms = numberOfAtoms;
	}

	public List<Integer> getAtomicNumbers() {
		return atomicNumbers;
	}

	public void setAtomicNumbers(List<Integer> atomicNumbers) {
		this.atomicNumbers = atomicNumbers;
	}

	public List<Vector3D> getThreeDCoords() {
		return threeDCoords;
	}

	public void setThreeDCoords(List<Vector3D> threeDCoords) {
		this.threeDCoords = threeDCoords;
	}

	public Double getEnergy() {
		return energy;
	}

	public void setEnergy(Double energy) {
		this.energy = energy;
	}

	public Double getMolecularMass() {
		return molecularMass;
	}

	public void setMolecularMass(Double molecularMass) {
		this.molecularMass = molecularMass;
	}

	public List<Double> getFrequencies() {
		return frequencies;
	}

	public void setFrequencies(List<Double> frequencies) {
		this.frequencies = frequencies;
	}

	public Integer getGroundStateDegeneracy() {
		return groundStateDegeneracy;
	}

	public void setGroundStateDegeneracy(Integer groundStateDegeneracy) {
		this.groundStateDegeneracy = groundStateDegeneracy;
	}

	public Double[] getRotationalConstants() {
		return rotationalConstants;
	}

	public void setRotationalConstants(Double[] rotationalConstants) {
		this.rotationalConstants = rotationalConstants;
	}

	public Double getStericEnergy() {
		return stericEnergy;
	}

	public void setStericEnergy(Double stericEnergy) {
		this.stericEnergy = stericEnergy;
	}
}
