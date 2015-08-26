package FastDAS;

import java.util.regex.Pattern;

/**
 * Created with IntelliJ IDEA.
 * User: uqlfearn
 * Date: 15/08/13
 * Time: 4:22 PM
 * To change this template use File | Settings | File Templates.
 */
public class KeyParam {

	public static final Object[] objectiveOptions = {"Inputs", "Outputs", "Reactions"};
	public static final Integer OBJECTIVE_INPUTS = 0;
	public static final Integer OBJECTIVE_OUTPUTS = 1;
	public static final Integer OBJECTIVE_REACTIONS = 2;

	public static final Object[] modeOptions = {"Minimise", "Maximise"};
	public static final Integer MODE_MINIMISE = 0;
	public static final Integer MODE_MAXIMISE = 1;

	public static final Object[] saveOptions = {"Full Cytoscape Network (SIF format, attribute files)",
			"Differential Cytoscape Network (Individual comparison networks, attribute files)", "Differential " +
			"Cytoscape Network (Union of comparison networks, attribute files)"};
	public static final Integer SAVE_FULL = 0;
	public static final Integer SAVE_DIFFERENTIAL_INDIVIDUAL = 1;
	public static final Integer SAVE_DIFFERENTIAL_UNION = 2;

	public static final int FREE = 0;
	public static final int INACTIVE = -1;
	public static final int ACTIVE = 1;

	public static final int CONTROL_CONDITION_INDEX = 0;

	public static final Pattern EXPERIMENT_NUMBER_PATTERN = Pattern.compile("EXP_(\\d+)");

	//UniProt Accession Patterns:
	//Pattern 1 : [A-N,R-Z][0-9][A-Z][A-Z, 0-9][A-Z, 0-9][0-9]
	public static final Pattern UNIPROT_PATTERN_ONE = Pattern.compile("([A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z, 0-9][0-9])");
	//Pattern 2 : [O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9]
	public static final Pattern UNIPROT_PATTERN_TWO = Pattern.compile("([O,P,Q][0-9][A-Z,0-9][A-Z,0-9][A-Z," +
			"0-9][0-9])");

	public static final Integer SETTER_USER_SPEC = 0;
	public static final Integer SETTER_ESET_FILL = 1;
	public static final Integer SETTER_COMPLEX_FILL = 2;

	public static final int PROBLEM_TYPE_SINGLE = 0;
	public static final int PROBLEM_TYPE_MULTI_TARGET = 1;

	public static final int TYPE_SINGLE_PROBLEM = 0;
	public static final int TYPE_TREESEARCH_PROBLEM = 1;

	public static final double SOLVECODE_DIAG_INFEASIBLE = -1.0;
	public static final double SOLVECODE_DIAG_OTHER = -2.0;

	public static final String CYTOSCAPE_FALSE = "-1";
	public static final String CYTOSCAPE_TRUE = "1";
	public static final String REACT_URL = "http://www.reactome.org/cgi-bin/eventbrowser_st_id?ST_ID=";

	public static final double SOLVER_INACTIVE = 0.0;
	public static final String VERSION_NUMBER = "0.1 ALPHA";
	public static final String REACTOME_ENTITYSET_FLAG = "Converted from EntitySet in Reactome";
	public static final String REACTOME_COMPLEX_ENTITYSET_CONTAINING_FLAG = "This Complex contains a Reactome " +
			"EntitySet as a component";

	public static final int OPEN = 1;
	public static final int CLOSED = 2;

	public static final int SPECIES_TYPE_IDX = 0;
	public static final int SPECIES_RDF_IDX = 1;
	public static final int SPECIES_NAME_IDX = 2;
	public static final int SPECIES_LOCATION_IDX = 3;
	public static final int SPECIES_PARTNUM_IDX = 4;
	public static final int SPECIES_IO_IDX = 5;
	public static final int SPECIES_DELETION_FLAG = 6;
	public static final int SPECIES_CONDITION_IDX_START = 7;

	public static final String NEWLINE = System.getProperty("line.separator");
	public static final String PATH_SEPARATOR = System.getProperty("file.separator");

	//Traversal limit for low-confidence outputs.
	public static final int TRAVERSAL_LIMIT = 100;

	public static final int TRAVERSAL_MED_CONF = 6;

	public static final boolean COMPUTE_RESPONSIVENESS_ON = false;
}
