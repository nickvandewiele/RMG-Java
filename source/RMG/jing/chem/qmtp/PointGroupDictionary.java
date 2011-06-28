package jing.chem.qmtp;

import java.util.HashMap;
import java.util.Map;

/**
 * Dictionary type with all the  symmetry groups supported in RMG
 * @author nmvdewie
 *
 */
public class PointGroupDictionary {

	//Map with the symmetry groups and their associated symmmetry numbers
	Map<String,Integer> library;

	//Map with the symmetry groups and a flag whether it is chiral or not
	Map<String,Boolean> chiralLibrary;



	private static PointGroupDictionary INSTANCE;

	private PointGroupDictionary(){
		library = new HashMap<String, Integer>();
		chiralLibrary = new HashMap<String, Boolean>(); 
		populateGroups();

		populateChiralityFlags();

	}

	private void populateChiralityFlags() {
		chiralLibrary.put("C1", true);
		chiralLibrary.put("Cs", false);
		chiralLibrary.put("Ci", false);
		chiralLibrary.put("C2", true);
		chiralLibrary.put("C3", true);
		chiralLibrary.put("C4",true);
		chiralLibrary.put("C5", true);
		chiralLibrary.put("C6",true);
		chiralLibrary.put("C7",true);
		chiralLibrary.put("C8",true);

		chiralLibrary.put("D2", true);
		chiralLibrary.put("D3", true);
		chiralLibrary.put("D4", true);
		chiralLibrary.put("D5", true);
		chiralLibrary.put("D6", true);
		chiralLibrary.put("D7", true);
		chiralLibrary.put("D8", true);

		chiralLibrary.put("C2v",false);
		chiralLibrary.put("C3v",false);
		chiralLibrary.put("C4v",false);
		chiralLibrary.put("C5v",false);
		chiralLibrary.put("C6v", false);
		chiralLibrary.put("C7v", false);
		chiralLibrary.put("C8v", false);

		chiralLibrary.put("C2h", false);
		chiralLibrary.put("C3h", false);
		chiralLibrary.put("C4h",false );
		chiralLibrary.put("C5h",false);
		chiralLibrary.put("C6h", false);
		chiralLibrary.put("C8h",false );

		chiralLibrary.put("D2h", false);
		chiralLibrary.put("D3h", false);
		chiralLibrary.put("D4h", false);
		chiralLibrary.put("D5h", false);
		chiralLibrary.put("D6h",false );
		chiralLibrary.put("D7h", false);
		chiralLibrary.put("D8h", false);

		chiralLibrary.put("D2d", false);
		chiralLibrary.put("D3d",false );
		chiralLibrary.put("D4d", false);
		chiralLibrary.put("D5d",false);
		chiralLibrary.put("D6d", false);
		chiralLibrary.put("D7d", false);
		chiralLibrary.put("D8d",false);

		chiralLibrary.put("S4", true);
		chiralLibrary.put("S6", true);
		chiralLibrary.put("S8", true);

		chiralLibrary.put("T",true );
		chiralLibrary.put("Th",false );
		chiralLibrary.put("Td",false);

		chiralLibrary.put("O", true);
		chiralLibrary.put("Oh",false );

		chiralLibrary.put("Cinfv",false );
		chiralLibrary.put("Dinfh",false );
		chiralLibrary.put("I", true);
		chiralLibrary.put("Ih",false );
		chiralLibrary.put("Kh",false );	

	}

	private void populateGroups() {
		library.put("C1", new Integer(1));
		library.put("Cs", new Integer(1));
		library.put("Ci", new Integer(1));
		library.put("C2", new Integer(2));
		library.put("C3", new Integer(3));
		library.put("C4", new Integer(4));
		library.put("C5", new Integer(5));
		library.put("C6", new Integer(6));
		library.put("C7", new Integer(7));
		library.put("C8", new Integer(8));

		library.put("D2", new Integer(4));
		library.put("D3", new Integer(6));
		library.put("D4", new Integer(8));
		library.put("D5", new Integer(10));
		library.put("D6", new Integer(12));
		library.put("D7", new Integer(14));
		library.put("D8", new Integer(16));

		library.put("C2v", new Integer(2));
		library.put("C3v", new Integer(3));
		library.put("C4v", new Integer(4));
		library.put("C5v", new Integer(5));
		library.put("C6v", new Integer(6));
		library.put("C7v", new Integer(7));
		library.put("C8v", new Integer(8));

		library.put("C2h", new Integer(2));
		library.put("C3h", new Integer(3));
		library.put("C4h", new Integer(4));
		library.put("C5h", new Integer(5));
		library.put("C6h", new Integer(6));
		library.put("C8h", new Integer(8));

		library.put("D2h", new Integer(4));
		library.put("D3h", new Integer(6));
		library.put("D4h", new Integer(8));
		library.put("D5h", new Integer(10));
		library.put("D6h", new Integer(12));
		library.put("D7h", new Integer(14));
		library.put("D8h", new Integer(16));

		library.put("D2d", new Integer(4));
		library.put("D3d", new Integer(6));
		library.put("D4d", new Integer(8));
		library.put("D5d", new Integer(10));
		library.put("D6d", new Integer(12));
		library.put("D7d", new Integer(14));
		library.put("D8d", new Integer(16));

		library.put("S4", new Integer(2));
		library.put("S6", new Integer(3));
		library.put("S8", new Integer(4));

		library.put("T", new Integer(12));
		library.put("Th", new Integer(12));
		library.put("Td", new Integer(12));

		library.put("O", new Integer(24));
		library.put("Oh", new Integer(24));

		library.put("Cinfv", new Integer(1));
		library.put("Dinfh", new Integer(2));
		library.put("I", new Integer(60));
		library.put("Ih", new Integer(60));
		library.put("Kh", new Integer(1));	

	}

	public boolean contains(String pointGroup){
		return library.keySet().contains(pointGroup);
	}
	public Integer get(String pointGroup){

		Integer symmNumber = library.get(pointGroup);

		if(symmNumber != null) return symmNumber;

		else return null;
	}

	public static synchronized PointGroupDictionary getInstance() {
		if(INSTANCE == null){
			INSTANCE = new PointGroupDictionary();
		}
		return INSTANCE;
	}

	public Boolean isChiral(String pointGroup) {
		Boolean chiral = chiralLibrary.get(pointGroup);

		if(chiral != null) return chiral;

		else return null;
	}
}
