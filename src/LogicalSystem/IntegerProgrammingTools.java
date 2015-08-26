package LogicalSystem;

import FastDAS.KeyParam;
import org.apache.log4j.Logger;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.Entity;
import org.biopax.paxtools.model.level3.Interaction;
import org.gnu.glpk.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created with IntelliJ IDEA.
 * User: uqlfearn
 * Date: 20/08/13
 * Time: 11:44 AM
 * To change this template use File | Settings | File Templates.
 */
public class IntegerProgrammingTools {

	private static org.apache.log4j.Logger log = Logger.getLogger(IntegerProgrammingTools.class);

	public static glp_prob getProblemLpSolveInstance(ProblemInstance instance) {
		try {
			glp_prob lp = null;
			//LpSolve lpInstance = null;
			HashMap<Integer, Object> integerObjectHashMap = instance.getIntegerObjectHashMap();
			HashMap<Object, Integer> objectIntegerHashMap = instance.getObjectIntegerHashMap();
			AdjList networkAdjList = instance.getNetworkAdjList();
			Model model = instance.getModel();
			//Create a new LpSolve instance. Number of rows (conditions) = 0 as we need to read the model into the
			//instance. Number of columns (variables) = number of entities in the hashmaps.
			lp = GLPK.glp_create_prob();
			//lpInstance = LpSolve.makeLp(0, integerObjectHashMap.size() + 1);
			//Set the solver verbosity - LpSolve.NORMAL is default. Consider setting to LpSolve.IMPORTANT (warnings
			//+ errors) or LpSolve.SEVERE (errors only) later.
			//lpInstance.setVerbose(LpSolve.SEVERE);
			//We're going to instantiate each variable in the hashmap one by one.
			GLPK.glp_add_cols(lp, integerObjectHashMap.keySet().size());
			for (Integer i : integerObjectHashMap.keySet()) {
				String name = "x"+i;
				//Set the name for the column.
				GLPK.glp_set_col_name(lp, i, ((Entity)integerObjectHashMap.get(i)).getRDFId());
				GLPK.glp_set_col_kind(lp, i, GLPK.GLP_BV);
			}
			//Set the lpInstance's behaviour for adding rows to expansion not overwrites.
			//lpInstance.setAddRowmode(true);
			//Now we populate the IP solver with the constraints responsible for activating each individual variable in
			//the problem.
			ArrayList<Integer> ilist = new ArrayList<Integer>();
			ArrayList<Integer> jlist = new ArrayList<Integer>();
			ArrayList<Double> vlist = new ArrayList<Double>();
			ilist.add(0);
			jlist.add(0);
			vlist.add(0.0);
			for (Object object : objectIntegerHashMap.keySet()) {
				//Interaction cases:
				if (object instanceof Interaction) {
					//Cast to Interaction
					Interaction currentInteraction = (Interaction) object;
					//Get the index of the interaction.
					int interactionIndex = objectIntegerHashMap.get(currentInteraction);
					//Get the left side of the reaction - these are the parent nodes of the interaction in the
					//network graph.
					HashSet<Integer> substrateSet = networkAdjList.getParents(interactionIndex);
					//Set up a new array to hold the constraints for the LHS of the reaction.
					double[] mainConstraint = new double[objectIntegerHashMap.size() + 1];
					//Case 1: No LHS. This -should- be caught before this step, but there's a warn here, just in
					//case...
					if (substrateSet.isEmpty()) {
						log.warn("WARN : Interaction " + currentInteraction.getRDFId().replaceAll(model.getXmlBase(),
								"")
								+ " with no left hand side in solver instantiation");
						//Case 2: LHS has one element.
					} else if (substrateSet.size() == 1) {
						//Constraints are of the form : x + y >= <value>. inequalityValue holds the value on the RHS.
						double inequalityValue = 0.0;
						int reactantIndex = substrateSet.iterator().next();
						//Standard reactant required for reaction to take place:
						if (reactantIndex >= 0) {
							mainConstraint[reactantIndex] = -1.0;
							//Or inhibiting interaction:
						} else {
							mainConstraint[Math.abs(reactantIndex)] = 1.0;
							inequalityValue = 1.0;
						}
						mainConstraint[interactionIndex] = 1.0;
						//Add value in.
						//Create new row:
						int rownum = GLPK.glp_add_rows(lp, 1);
						GLPK.glp_set_row_bnds(lp, rownum, GLPKConstants.GLP_FX, inequalityValue, inequalityValue);
						for (int i = 0; i < mainConstraint.length;i++) {
							if (mainConstraint[i] != 0.0) {
								ilist.add(rownum);
								jlist.add(i);
								vlist.add(mainConstraint[i]);
							}
						}
						//lpInstance.addConstraint(mainConstraint, LpSolve.EQ,
						//		inequalityValue);
						//Case 3: LHS > 1 element.
					} else {
						double inequalityValue = 1.0;
						for (int substrateIndex : substrateSet) {
							double[] elementConstraint = new double[objectIntegerHashMap.size() + 1];
							//Standard reactant required for reaction to take place:
							if (substrateIndex >= 0) {
								//lpInstance.addConstraint(elementConstraint, LpSolve.GE, 0.0);
								//elementConstraint[substrateIndex] = 1.0;
								//elementConstraint[interactionIndex] = -1.0;
								int rownum = GLPK.glp_add_rows(lp, 1);
								ilist.add(rownum);
								jlist.add(substrateIndex);
								vlist.add(1.0);
								ilist.add(rownum);
								jlist.add(interactionIndex);
								vlist.add(-1.0);
								GLPK.glp_set_row_bnds(lp, rownum, GLPKConstants.GLP_LO, 0.0, 0.0);
								mainConstraint[substrateIndex] = -1.0;
								inequalityValue = inequalityValue - 1.0;
								//Or inhibiting:
							} else {
								//elementConstraint[Math.abs(substrateIndex)] = -1.0;
								//elementConstraint[interactionIndex] = -1.0;
								//lpInstance.addConstraint(elementConstraint, LpSolve.GE, -1.0);
								int rownum = GLPK.glp_add_rows(lp,1);
								ilist.add(rownum);
								jlist.add(Math.abs(substrateIndex));
								vlist.add(-1.0);
								ilist.add(rownum);
								jlist.add(interactionIndex);
								vlist.add(-1.0);
								GLPK.glp_set_row_bnds(lp, rownum, GLPKConstants.GLP_LO, -1.0, -1.0);
								mainConstraint[Math.abs(substrateIndex)] = 1.0;
							}
						}
						mainConstraint[interactionIndex] = 1.0;
						//lpInstance.addConstraint(mainConstraint, LpSolve.GE, inequalityValue);
						int rownum = GLPK.glp_add_rows(lp,1);
						for (int i = 0; i < mainConstraint.length;i++) {
							if (mainConstraint[i] != 0.0) {
								ilist.add(rownum);
								jlist.add(i);
								vlist.add(mainConstraint[i]);
							}
						}
						GLPK.glp_set_row_bnds(lp, rownum, GLPKConstants.GLP_LO, inequalityValue, inequalityValue);
					}
					//Else element is a physical entity participating in the system:
				} else {
					int speciesIndex = objectIntegerHashMap.get(object);
					//Get the left side of the reaction - these are the parent nodes of the interaction in the
					//network graph.
					HashSet<Integer> activatingEvents = networkAdjList.getParents(speciesIndex);
					//Set up a new array to hold the constraints for the LHS of the reaction.
					double[] mainConstraint = new double[objectIntegerHashMap.size() + 1];
					//Case 1: No activating events. This is ok - this just means that the species is an input into the
					//network, i.e., nothing produces it.
					if (activatingEvents.isEmpty()) {
						//Skip forward, species is an input.
						continue;
						//Case 2: species is produced by one reaction.
					} else if (activatingEvents.size() == 1) {
						//Constraints are of the form : x + y >= <value>. inequalityValue holds the value on the RHS.
						double inequalityValue = 0.0;
						int reactionIndex = activatingEvents.iterator().next();
						//Standard reactant required for reaction to take place:
						if (reactionIndex >= 0) {
							mainConstraint[reactionIndex] = 1.0;
							//Reactions should not be inhibiting species (we use a mechanistic representation of
							// inhibition
							//due to a lack of relative reaction kinetics that would allow us to reasonably discuss
							//competitive 'inhibition' (e.g. rapid degradation halting signalling). Biologically this
							// happens,
							//but the data is not available in a usable form.
						} else {
							log.debug("Reaction " +
									((Entity) integerObjectHashMap.get(reactionIndex)).getRDFId().replaceAll(model
											.getXmlBase(), "")
									+ " inhibiting species " +
									((Entity) integerObjectHashMap.get(speciesIndex)).getRDFId().replaceAll(model
											.getXmlBase(), ""));
						}
						mainConstraint[speciesIndex] = -1.0;
						//Add the value in.
						//lpInstance.addConstraint(mainConstraint, LpSolve.EQ, 0.0);
						int rownum = GLPK.glp_add_rows(lp,1);
						for (int i = 0; i < mainConstraint.length;i++) {
							if (mainConstraint[i] != 0.0) {
								ilist.add(rownum);
								jlist.add(i);
								vlist.add(mainConstraint[i]);
							}
						}
						GLPK.glp_set_row_bnds(lp, rownum, GLPKConstants.GLP_FX, 0.0, 0.0);
						//Case 3: LHS > 1 element.
					} else {
						for (int reactionIndex : activatingEvents) {
							double[] elementConstraint = new double[objectIntegerHashMap.size() + 1];
							//Standard reactant required for reaction to take place:
							if (reactionIndex >= 0) {
								//elementConstraint[reactionIndex] = -1.0;
								//elementConstraint[speciesIndex] = 1.0;
								//lpInstance.addConstraint(elementConstraint, LpSolve.GE, 0.0);
								int rownum = GLPK.glp_add_rows(lp, 1);
								ilist.add(rownum);
								jlist.add(reactionIndex);
								vlist.add(-1.0);
								ilist.add(rownum);
								jlist.add(speciesIndex);
								vlist.add(1.0);
								GLPK.glp_set_row_bnds(lp, rownum, GLPKConstants.GLP_LO, 0.0, 0.0);
								mainConstraint[reactionIndex] = 1.0;
								//Reactions should not be inhibiting species (we use a mechanistic representation of
								// inhibition
								//due to a lack of relative reaction kinetics that would allow us to reasonably discuss
								//competitive 'inhibition' (e.g. rapid degradation halting signalling). Biologically
								// this happens,
								//but the data is not available in a usable form.
							} else {
								log.debug("Reaction " +
										((Entity) integerObjectHashMap.get(reactionIndex)).getRDFId().replaceAll(model
												.getXmlBase(), "")
										+ " inhibiting species " +
										((Entity) integerObjectHashMap.get(speciesIndex)).getRDFId().replaceAll(model
												.getXmlBase(), ""));
							}
						}
						mainConstraint[speciesIndex] = -1.0;
						//Reactions can't inhibit species, so the RHS of the inequality is 0.0 (there's no integer
						// terms
						//on the LHS of the inequality to be shifted over during instantiation).
						//lpInstance.addConstraint(mainConstraint, LpSolve.GE, 0.0);
						int rownum = GLPK.glp_add_rows(lp,1);
						for (int i = 0; i < mainConstraint.length;i++) {
							if (mainConstraint[i] != 0.0) {
								ilist.add(rownum);
								jlist.add(i);
								vlist.add(mainConstraint[i]);
							}
						}
						GLPK.glp_set_row_bnds(lp, rownum, GLPKConstants.GLP_LO, 0.0, 0.0);
					}
				}
			}
			SWIGTYPE_p_int iarray = GLPK.new_intArray(ilist.size());
			SWIGTYPE_p_int jarray = GLPK.new_intArray(jlist.size());
			SWIGTYPE_p_double varray = GLPK.new_doubleArray(vlist.size());
			for (int i = 1; i < ilist.size(); i++) {
				GLPK.intArray_setitem(iarray, i, ilist.get(i));
				GLPK.intArray_setitem(jarray, i, jlist.get(i));
				GLPK.doubleArray_setitem(varray, i, vlist.get(i));
			}
			GLPK.glp_load_matrix(lp, ilist.size()-1, iarray, jarray, varray);
			GLPK.delete_doubleArray(varray);
			GLPK.delete_intArray(iarray);
			GLPK.delete_intArray(jarray);
			return lp;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	public static glp_prob setObjectiveFunction(glp_prob originalLP, HashSet<Integer> objective) {
		glp_prob newLP = null;
		try {
			newLP = GLPK.glp_create_prob();
			GLPK.glp_copy_prob(newLP, originalLP, GLPKConstants.GLP_ON);
			//derivativeInstance = originalLP.copyLp();
			//double[] objectiveArray = new double[originalLP.getNcolumns() + 1];
			for (Integer i : objective) {
				if (i < 0) {
					//TODO:Check this....
					//objectiveArray[Math.abs(i)] = -1.0;
					GLPK.glp_set_obj_coef(newLP, Math.abs(i), -1.0);
				} else {
					GLPK.glp_set_obj_coef(newLP, Math.abs(i), 1.0);
					//objectiveArray[Math.abs(i)] = 1.0;
				}
			}
			//derivativeInstance.setObjFn(objectiveArray);
		} catch (Exception e) {
			log.error("Error in setting objective.", e);
		} finally {
			return newLP;
		}
	}

	public static ArrayList<glp_prob> setConstraints(glp_prob mainLP, HashMap<Object,
			ArrayList<Integer>> instanceSettingsMap, int conditionCount, HashMap<Object, Integer>
															objectIntegerHashMap) {
		ArrayList<glp_prob> lpSolves = new ArrayList<glp_prob>(conditionCount);
		try {
			for (int i = 0; i < conditionCount; i++) {
				glp_prob currentLP = GLPK.glp_create_prob();
				GLPK.glp_copy_prob(currentLP, mainLP, GLPK.GLP_ON);
				for (Object obj : objectIntegerHashMap.keySet()) {
					if (!instanceSettingsMap.containsKey(obj)) {
						continue;
					}
					Integer index = objectIntegerHashMap.get(obj);
					if (index==null) {System.out.println(((Entity)obj).getRDFId() + " FAIL HERE"); System.out.println
							(index); System
							.exit(-1);}
					if (instanceSettingsMap.get(obj).get(i).equals(KeyParam.ACTIVE)) {
						//Set the constraint for that variable to be active.
						GLPK.glp_set_col_bnds(currentLP, index, GLPK.GLP_FX, 1.0, 1.0);
						//currentLP.setBounds(index, 1.0, 1.0);
					} else if (instanceSettingsMap.get(obj).get(i).equals(KeyParam.INACTIVE)) {
						//Set the constraint for that variable to be inactive.
						GLPK.glp_set_col_bnds(currentLP, index, GLPK.GLP_FX, 0.0, 0.0);
						//currentLP.setBounds(index, 0.0, 0.0);
					} else {
						//Free, don't worry about it.
					}
				}
				lpSolves.add(currentLP);
			}
		} catch (Exception e) {
			log.error("ERROR: Failed generation of condition LP instances", e);
		} finally {
			return lpSolves;
		}
	}

	public static void cut(glp_prob lpInstance, HashSet<Integer> newCutSet) {
		for (Integer index : newCutSet) {
			try {
				GLPK.glp_set_col_bnds(lpInstance, index, GLPKConstants.GLP_FX, 0.0, 0.0);
				//lpInstance.setBounds(index, 0.0, 0.0);
			} catch (Exception e) {
				log.error("ERROR: Failed to set cut condition", e);
			}
		}
	}
}
