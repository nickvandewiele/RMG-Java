package jing.chem.qmtp;

public class PointGroup {
	
	String pointGroup;
	
	//associated Symmetry Number
	Integer symmetryNumber;
	
	//chiral flag
	Boolean chiral;
	
	//linear flag
	Boolean linear;
	
	public PointGroup(String pointGroup){
		this.pointGroup = pointGroup;
		
		this.symmetryNumber = PointGroupDictionary.getInstance().get(pointGroup);
		this.chiral = PointGroupDictionary.getInstance().isChiral(pointGroup);
		
		//determine linearity from 3D-geometry; changed to correctly consider linear ketene radical case
		if (pointGroup.equals("Cinfv")||pointGroup.equals("Dinfh")) linear = true;
		else linear = false;
	}
	
	@Override
	public String toString(){
		return pointGroup;
	}
	
	public boolean equals(String group){
		return pointGroup.equals(group);
	}
	
	public Boolean isLinear(){
		return linear;
	}
}
