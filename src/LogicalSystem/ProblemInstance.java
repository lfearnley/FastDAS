package LogicalSystem;

import FastDAS.KeyParam;
import IO.TextOutputHandler;
import Testing.TestHarness;
import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.security.Key;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.*;
import java.util.regex.Matcher;
import org.gnu.glpk.*;

/**
 * Created with IntelliJ IDEA.
 * User: uqlfearn
 * Date: 16/08/13
 * Time: 3:39 PM
 * To change this template use File | Settings | File Templates.
 */
public class ProblemInstance {

	private org.apache.log4j.Logger log = Logger.getLogger(ProblemInstance.class);

	private Model model;
	private HashMap<Integer, Object> integerObjectHashMap;
	private HashMap<Object, Integer> objectIntegerHashMap;
	private AdjList networkAdjList = new AdjList();
	private HashMap<Object, ArrayList<Integer>> instanceSettingsMap;
	private HashMap<Object, Integer> setterMap;
	private int conditionCount = 1;
	private int type = KeyParam.PROBLEM_TYPE_SINGLE;
	private int mode = KeyParam.TYPE_SINGLE_PROBLEM;
	private int optMode = KeyParam.MODE_MAXIMISE;
	private int outputMode = KeyParam.SAVE_DIFFERENTIAL_INDIVIDUAL;
	private HashSet<Integer> objectiveSet;
	private HashSet<Integer> reactionSet;
	private ArrayList<double[]> solutions;
	private File outputDirectory;
	private boolean saveMidConf;
	private boolean saveLowConf;
	private String descriptor = "";
	private glp_prob mainLP;


	private LinkedBlockingQueue<ArrayList<Integer>> solsizes;

	private boolean enumerateStats = false;

	public ProblemInstance(Model model, Integer objective, Integer optMode, Integer outputMode, File outputDirectory,
						   boolean saveMidConf, boolean saveLowConf, String descriptor) {
		this.model = model;
		this.conditionCount = getConditionCount();
		this.integerObjectHashMap = new HashMap<Integer, Object>();
		this.objectIntegerHashMap = new HashMap<Object, Integer>();
		this.optMode = optMode;
		this.saveMidConf = saveMidConf;
		this.saveLowConf = saveLowConf;
		this.descriptor = descriptor;
		System.out.println("INTO INSTANCE POP");
		populateProblemInstance();
		System.out.println("SETTINGS MAP POP");
		if (!KeyParam.COMPUTE_RESPONSIVENESS_ON) {
			populateInstanceSettingsMap();
			System.out.println("THROUGH INSTANCE SETTINGS");
			this.reactionSet = new HashSet<Integer>();
			for (Object obj : objectIntegerHashMap.keySet()) {
				if (obj instanceof Conversion) {
					reactionSet.add(objectIntegerHashMap.get(obj));
					if (!obj.equals(integerObjectHashMap.get(objectIntegerHashMap.get(obj)))) {
						System.err.println("HASHMAP DOESN'T MATCH");
						System.exit(-1);
					}
				}
			}
			System.out.println("THROUGH OIMAP");
			this.outputDirectory = outputDirectory;
			this.outputMode = outputMode;
			if (objective.equals(KeyParam.OBJECTIVE_INPUTS)) {
				objectiveSet = networkAdjList.getInputs();
			} else if (objective.equals(KeyParam.OBJECTIVE_OUTPUTS)) {
				objectiveSet = networkAdjList.getOutputs();
			} else if (objective.equals(KeyParam.OBJECTIVE_REACTIONS)) {
				System.out.println("SETTING REACTIONS");
				objectiveSet = new HashSet<Integer>(reactionSet);
			} else {
				System.out.println("COBJ EXCEPTION");
				throw new UnsupportedOperationException("CUSTOM OBJECTIVES NOT WRITTEN YET.");
			}
			try {
				System.out.println("XSPINNING UP");
				populateMainLP();
				System.out.println("SOLVERSPIN");
				spinUpSolvers();
			} catch (Exception e) {
				log.error("ERROR: solver error", e);
				System.exit(-1);
			}
		} else {
			this.outputDirectory = outputDirectory;
			BufferedWriter writer = null;
			BufferedWriter failwriter = null;
			try {
				solsizes = new LinkedBlockingQueue<ArrayList<Integer>>();
				writer = new BufferedWriter(new FileWriter(new File(this.getOutputDirectory() +
						KeyParam.PATH_SEPARATOR + "output.txt")));
				writer.write("HIGH\tMED\tLOW"+KeyParam.NEWLINE);
				this.conditionCount = 2;
				this.reactionSet = new HashSet<Integer>();
				for (Object obj : objectIntegerHashMap.keySet()) {
					if (obj instanceof Conversion) {
						reactionSet.add(objectIntegerHashMap.get(obj));
						if (!obj.equals(integerObjectHashMap.get(objectIntegerHashMap.get(obj)))) {
							System.err.println("HASHMAP DOESN'T MATCH");
							System.exit(-1);
						}
					}
				}
				System.out.println("THROUGH OIMAP");
				this.outputDirectory = outputDirectory;
				this.outputMode = outputMode;
				if (objective.equals(KeyParam.OBJECTIVE_INPUTS)) {
					objectiveSet = networkAdjList.getInputs();
				} else if (objective.equals(KeyParam.OBJECTIVE_OUTPUTS)) {
					objectiveSet = networkAdjList.getOutputs();
				} else if (objective.equals(KeyParam.OBJECTIVE_REACTIONS)) {
					System.out.println("SETTING REACTIONS");
					objectiveSet = new HashSet<Integer>(reactionSet);
				} else {
					System.out.println("COBJ EXCEPTION");
					throw new UnsupportedOperationException("CUSTOM OBJECTIVES NOT WRITTEN YET.");
				}
				int inputcount = this.getNetworkAdjList().getInputs().size();
				int counter = 0;
				int runcount = 0;
				ArrayList<String> failedlist = new ArrayList<String>();
				for (int i = 0; i < 14; i++) {
					failedlist.add("");
				}
				this.instanceSettingsMap = new HashMap<Object, ArrayList<Integer>>();
				for (Object obj : objectIntegerHashMap.keySet()) {
					if (!instanceSettingsMap.containsKey(obj)) {
						ArrayList<Integer> integerList = new ArrayList<Integer>(conditionCount);
						for (int i = 0; i < conditionCount; i++) {
							integerList.add(KeyParam.FREE);
						}
						instanceSettingsMap.put(obj, integerList);
					}
				}
				this.setterMap = new HashMap<Object, Integer>();
				for (Integer input : this.getNetworkAdjList().getInputs()) {
					if (this.getNetworkAdjList().isOutput(input)) {
						continue;
					}
					runcount++;
					System.out.print(counter);
					System.out.print("\\");
					System.out.println(inputcount);
					counter++;
					Object inputObject = integerObjectHashMap.get(Math.abs(input));
					if (inputObject == null) {
						System.out.println("GOT NULL OBJECT");
						System.out.println(input);
					}
					System.out.println(((Entity) inputObject).getRDFId());
					instanceSettingsMap.get(inputObject).set(0, KeyParam.INACTIVE);
					instanceSettingsMap.get(inputObject).set(1, KeyParam.ACTIVE);
					try {
						System.out.println(new SimpleDateFormat("hhmm").format(new Date()));
						populateMainLP();
						//System.out.println("SPINNING UP");
						spinUpSolvers();
					} catch (Exception e) {
						log.error("ERROR: solver error", e);
						System.exit(-1);
					}
					StringBuilder output = new StringBuilder();
					int errorcode = 0;
					ArrayList<Integer> sizes = solsizes.poll();
					for (int i = 0; i < 3; i++) {
						if (sizes.get(i) < 0) {
							errorcode = Math.abs(sizes.get(i));
						}
						output.append(sizes.get(i).toString());
						output.append("\t");
					}
					if (errorcode == 0) {
						writer.write(output + KeyParam.NEWLINE);
					} else {
						String ecString = failedlist.get(errorcode) + ((Entity) inputObject).getRDFId() + KeyParam.NEWLINE;

						failedlist.set(errorcode, ecString);
					}
					instanceSettingsMap.get(inputObject).set(0, KeyParam.FREE);
					instanceSettingsMap.get(inputObject).set(1, KeyParam.FREE);
					if (runcount % 20 == 0) {
						System.gc();
						writer.flush();
					}
				}
			} catch (Exception e) {
				e.printStackTrace();
			} finally {
				IOUtils.closeQuietly(writer);
				IOUtils.closeQuietly(failwriter);
			}
		}
	}

	private int getConditionCount() {
		PhysicalEntity conditionCountObject = (PhysicalEntity) model.getByID("ConditionCount");
		if (conditionCountObject==null) {
			System.out.println("NULL OBJECT");
			return 1;
		}
		for (String s : conditionCountObject.getComment()) {
			if (!Character.isDigit(s.charAt(0))) {
				continue;
			} else {
				System.out.println("GOT CONDITION COUNT OF " + s);
				return Integer.parseInt(s);
			}
		}
		System.out.println("RETURNING 1");
		return 1;
	}

	private void populateInstanceSettingsMap() {
		this.setterMap = new HashMap<Object, Integer>();
		this.instanceSettingsMap = new HashMap<Object, ArrayList<Integer>>();
		System.out.println(conditionCount);
		for (Object obj : objectIntegerHashMap.keySet()) {
			if (!instanceSettingsMap.containsKey(obj)) {
				ArrayList<Integer> integerList = new ArrayList<Integer>(conditionCount);
				for (int i = 0; i < conditionCount; i++) {
					integerList.add(KeyParam.FREE);
				}
				instanceSettingsMap.put(obj, integerList);
			}
		}
		System.out.println("FILLED INSTANCESETTINGSMAP ConditionCount = " + conditionCount);
		for (Entity entity : model.getObjects(Entity.class)) {
			//Only process PhysicalEntities or Conversions.
			if (!(entity instanceof PhysicalEntity || entity instanceof Conversion)) {
				continue;
			}
			ArrayList<Integer> experimentalSettings = new ArrayList<Integer>(conditionCount);
			for (int i = 0; i < conditionCount; i++) {
				experimentalSettings.add(KeyParam.FREE);
			}
			//Flag to note that the entity being processed has experimental settings tagged in.
			boolean hasExperimentalSpec = false;
			//Run through each comment looking for the experimental flags.
			for (String commentString : entity.getComment()) {
				//Control/Condition '0'. This is the condition to which all comparisons are made.
				if (commentString.contains("CONTROL")) {
					if (commentString.contains("INACTIVE")) {
						hasExperimentalSpec = true;
						experimentalSettings.set(KeyParam.CONTROL_CONDITION_INDEX, KeyParam.INACTIVE);
						if (instanceSettingsMap.containsKey(entity)) {
							setterMap.put(entity, KeyParam.SETTER_USER_SPEC);
							instanceSettingsMap.put(entity, experimentalSettings);
						}
					} else if (commentString.contains("ACTIVE")) {
						hasExperimentalSpec = true;
						experimentalSettings.set(KeyParam.CONTROL_CONDITION_INDEX, KeyParam.ACTIVE);
						if (instanceSettingsMap.containsKey(entity)) {
							setterMap.put(entity, KeyParam.SETTER_USER_SPEC);
							instanceSettingsMap.put(entity, experimentalSettings);
						}
					} else {
						log.warn("Entity " + entity.getRDFId().replaceAll(model.getXmlBase(), "") + " failed control" +
								" " +
								"condition setting due to malformed specification.");
					}
					//Subsequent conditions. These are denoted using EXP_<conditionNumber> <setting>.
				} else if (commentString.contains("EXP_")) {
					Matcher matcher = KeyParam.EXPERIMENT_NUMBER_PATTERN.matcher(commentString);
					if (matcher.find()) {
						int val = Integer.parseInt(matcher.group(1));
						if (commentString.contains("INACTIVE")) {
							hasExperimentalSpec = true;
							experimentalSettings.set(val, KeyParam.INACTIVE);
							if (instanceSettingsMap.containsKey(entity)) {
								setterMap.put(entity, KeyParam.SETTER_USER_SPEC);
								instanceSettingsMap.put(entity, experimentalSettings);
							}
						} else if (commentString.contains("ACTIVE")) {
							hasExperimentalSpec = true;
							experimentalSettings.set(val, KeyParam.ACTIVE);
							if (instanceSettingsMap.containsKey(entity)) {
								setterMap.put(entity, KeyParam.SETTER_USER_SPEC);
								instanceSettingsMap.put(entity, experimentalSettings);
							}
						} else {
							log.warn("Entity " + entity.getRDFId().replaceAll(model.getXmlBase(), "") + " failed " +
									"experimental condition setting due to malformed specification of EXP_" + val);
						}
					} else {
						log.warn("Entity " + entity.getRDFId().replaceAll(model.getXmlBase(), "") + " failed " +
								"experimental " +
								"condition setting due to malformed specification.");
					}
				}
			}
			if (hasExperimentalSpec) {
				instanceSettingsMap.put(entity, experimentalSettings);
			}
			/*if (hasExperimentalSpec && entity instanceof PhysicalEntity) {
				PhysicalEntity pe = (PhysicalEntity) entity;
                //Bubble up to complexes and PhysicalEntities.
                LinkedBlockingQueue<Complex> complexQueue = new LinkedBlockingQueue<Complex>();
                LinkedBlockingQueue<PhysicalEntity> setQueue = new LinkedBlockingQueue<PhysicalEntity>();
                complexQueue.addAll(pe.getComponentOf());
                setQueue.addAll(pe.getMemberPhysicalEntityOf());
                HashSet<PhysicalEntity> doneset = new HashSet<PhysicalEntity>();
                while (!(complexQueue.isEmpty() && setQueue.isEmpty())) {
                    PhysicalEntity current;
                    if (complexQueue.isEmpty()) {
                        current = complexQueue.poll();
                    } else {
                        current = setQueue.poll();
                    }
                    if (current == null) {
                        continue;
                    }
                    if (objectIntegerHashMap.containsKey(current)) {
                        //If the node has been user-specified, terminate bubble at this point (we don't want to override
                        //user settings):
                        if (setterMap.containsKey(current) && setterMap.get(current).equals(KeyParam
                        .SETTER_USER_SPEC)) {
                            continue;
                        }
                        //If the node is an input, set it.
                        if (networkAdjList.isInput(objectIntegerHashMap.get(current))) {
                            instanceSettingsMap.put(current, experimentalSettings);
                            if (current instanceof Complex) {
                                setterMap.put(current, KeyParam.SETTER_COMPLEX_FILL);
                            } else {
                                setterMap.put(current, KeyParam.SETTER_ESET_FILL);
                            }
                        }
                    }
                    //Bubble up:
                    for (PhysicalEntity parent : current.getComponentOf()) {
                        if (!doneset.contains(parent)) {
                            complexQueue.add((Complex)parent);
                        }
                    }
                    for (PhysicalEntity set : current.getMemberPhysicalEntityOf()) {
                        if (!doneset.contains(set)) {
                            setQueue.add(set);
                        }
                    }
                    for (PhysicalEntity set : current.getMemberPhysicalEntity()) {
                        if (!doneset.contains(set)) {
                            setQueue.add(set);
                        }
                    }
                }
            }*/
		}
	}

	public HashMap<Integer, Object> getIntegerObjectHashMap() {
		return integerObjectHashMap;
	}

	public HashMap<Object, Integer> getObjectIntegerHashMap() {
		return objectIntegerHashMap;
	}

	public HashSet<Integer> getReactionSet() {
		return reactionSet;
	}

	public AdjList getNetworkAdjList() {
		return networkAdjList;
	}

	public Model getModel() {
		return model;
	}

	private void populateProblemInstance() {
		LinkedBlockingQueue<Control> controlQueue = new LinkedBlockingQueue<Control>();
		for (Interaction interaction : model.getObjects(Interaction.class)) {
			if (interaction instanceof Conversion) {
				Conversion conv = (Conversion) interaction;
				int convindex;
				if (objectIntegerHashMap.containsKey(conv)) {
					convindex = objectIntegerHashMap.get(conv);
				} else {
					convindex = addObjectToMaps(conv);
				}
				for (PhysicalEntity pe : conv.getLeft()) {
					if (!objectIntegerHashMap.containsKey(pe)) {
						addObjectToMaps((pe));
					}
					networkAdjList.addEdge(objectIntegerHashMap.get(pe), convindex);
				}
				for (PhysicalEntity pe : conv.getRight()) {
					if (!objectIntegerHashMap.containsKey(pe)) {
						addObjectToMaps(pe);
					}
					networkAdjList.addEdge(convindex, objectIntegerHashMap.get(pe));
				}
				for (Control c : conv.getControlledOf()) {
					ControlType controlType = c.getControlType();
					for (Controller controller : c.getController()) {
						if (controller instanceof PhysicalEntity) {
							PhysicalEntity pe = (PhysicalEntity) controller;
							if (!objectIntegerHashMap.containsKey(pe)) {
								addObjectToMaps(pe);
							}
							if (controlType.equals(ControlType.INHIBITION) || controlType.equals(ControlType
									.INHIBITION_OTHER)
									|| controlType.equals(ControlType.INHIBITION_ALLOSTERIC)
									|| controlType.equals(ControlType.INHIBITION_COMPETITIVE)
									|| controlType.equals(ControlType.INHIBITION_IRREVERSIBLE)
									|| controlType.equals(ControlType.INHIBITION_NONCOMPETITIVE)
									|| controlType.equals(ControlType.INHIBITION_UNCOMPETITIVE)
									|| controlType.equals(ControlType.INHIBITION_UNKMECH)) {
								networkAdjList.addEdge((-1 * objectIntegerHashMap.get(pe)), convindex);
							} else {
								networkAdjList.addEdge(objectIntegerHashMap.get(pe), convindex);
							}
						}
					}
				}
			} else {
				//log.warn("Not capturing " + interaction.getRDFId().replaceAll(model.getXmlBase(), ""));
			}
		}
	}

	private int addObjectToMaps(Object object) {
		//We want to index to IP solver column - this means count from 1, not 0, hence +1 offset on hashmap.
		objectIntegerHashMap.put(object, objectIntegerHashMap.size() + 1);
		integerObjectHashMap.put(objectIntegerHashMap.get(object), object);
		return objectIntegerHashMap.get(object);
	}

	private void populateMainLP() {
		this.mainLP = IntegerProgrammingTools.getProblemLpSolveInstance(this);
		//Set the objective function.
		mainLP = IntegerProgrammingTools.setObjectiveFunction(mainLP, objectiveSet);
		if (this.optMode == KeyParam.MODE_MAXIMISE) {
			//System.out.println("MAXMODE");
			GLPK.glp_set_obj_dir(mainLP, GLPKConstants.GLP_MAX);
		} else {
			GLPK.glp_set_obj_dir(mainLP, GLPKConstants.GLP_MIN);
		}
		//Set the constraints
		//mainLP.writeLp("main.lp");

	}

	private void spinUpSolvers() {
		//Get LP instance (constraints etc).
		/*for (int i = 0; i < lpArray.size(); i++) {
			LpSolve current = lpArray.get(i);
			current.writeLp("lp_" + i + ".lp");
		}*/
		ArrayList<glp_prob> lpArray = IntegerProgrammingTools.setConstraints(mainLP, instanceSettingsMap,
				conditionCount,
				objectIntegerHashMap);
		ExecutorService pool = null;
		try {
			//Create a cached thread pool for execution/dispatch:
			pool = Executors.newFixedThreadPool(1);
			//Single problem:
			if (this.mode == KeyParam.TYPE_SINGLE_PROBLEM) {
				//System.out.println("INTO MODE");
				//Working only on one instance of the network (as opposed to all outputs of the network individually):
				if (this.type == KeyParam.PROBLEM_TYPE_SINGLE) {
					ArrayList<Future<double[]>> results = new ArrayList<Future<double[]>>(lpArray.size());
					//Submit all callable LpSolves to the executor
					for (int i = 0; i < lpArray.size(); i++) {
						results.add(pool.submit(new LpSolveCallable(lpArray.get(i))));
					}
					pool.shutdown();
					if (!pool.awaitTermination(100, TimeUnit.SECONDS)) {
						pool.shutdownNow();
						if (!pool.awaitTermination(60, TimeUnit.SECONDS)) {
							System.err.println("FAILURE TO TERMINATE");
						}
					}
					//System.out.println("INTO OUTPUT");
					if (!KeyParam.COMPUTE_RESPONSIVENESS_ON) {
						//Write the results (spooling off the solver queue):
						if (this.outputMode == KeyParam.SAVE_DIFFERENTIAL_INDIVIDUAL) {
							IO.NetworkOutputHandler.writeDifferentialGraph(results, this, false);
						} else if (this.outputMode == KeyParam.SAVE_DIFFERENTIAL_UNION) {
							IO.NetworkOutputHandler.writeDifferentialGraph(results, this, true);
						} else if (this.outputMode == KeyParam.SAVE_FULL) {
							IO.NetworkOutputHandler.writeFullGraph(results, this);
						}
					} else {
						double[] control = results.get(0).get();
						double[] condition = results.get(1).get();
						if (control.length == 1) {
							System.err.println("CONTROL INFEASIBLE: " + control[0]);
						}
						if (condition.length == 1) {
							System.err.println("CONDITION INFEASIBLE: " + condition[0]);
						}
						if (condition.length != 1 && control.length != 1) {
							ArrayList<Integer> difference = Util.getDifference(control, condition,
									this);
							System.out.print(difference.get(0));
							System.out.print("\t");
							System.out.print(difference.get(1));
							System.out.print("\t");
							System.out.println(difference.get(2));
							solsizes.put(difference);
						} else if (control.length == 1) {
							ArrayList<Integer> al = new ArrayList<Integer>();
							al.add((new Double(control[0]).intValue()));
 							solsizes.put(al);
						} else if (condition.length == 1) {
							ArrayList<Integer> al = new ArrayList<Integer>();
							al.add((new Double(condition[0]).intValue()));
							solsizes.put(al);
						}
					}
				}
				//Branch-and-bound algorithm for integer cut enumeration
			} else if (this.mode == KeyParam.TYPE_TREESEARCH_PROBLEM) {
				if (this.type == KeyParam.PROBLEM_TYPE_SINGLE) {
					LinkedBlockingQueue<double[]> resultQueue = new LinkedBlockingQueue<double[]>();
					for (int i = 0; i < lpArray.size(); i++) {
						//Results are ordered by objective function size. Each objective function/optimal solution
						//count has an associated list of solutions, comprising a HashSet of integers describing
						//nodes active in the network, and a HashSet of integers describing nodes cut from the network.
						HashMap<Integer, ArrayList<Tuple2<HashSet<Integer>, HashSet<Integer>>>> resultMap
								= new HashMap<Integer, ArrayList<Tuple2<HashSet<Integer>, HashSet<Integer>>>>();
						//HashSet of HashSets - each component HashSet represents a set of integer cuts that have been
						//analysed/queued up for lpSolve instance i.
						HashSet<HashSet<Integer>> completedCutSets = new HashSet<HashSet<Integer>>();
						//Create and store the empty cutset.
						HashSet<Integer> parentCutSet = new HashSet<Integer>();
						completedCutSets.add(parentCutSet);
						LinkedBlockingQueue<Future<Tuple2<double[], HashSet<Integer>>>> callableQueue =
								new LinkedBlockingQueue<Future<Tuple2<double[], HashSet<Integer>>>>(lpArray.size());
						//Queue the first problem.
						glp_prob copyProb = GLPK.glp_create_prob();
						GLPK.glp_copy_prob(copyProb, lpArray.get(i), GLPKConstants.GLP_ON);
						callableQueue.add(pool.submit(new LpSolveCutSetCallable(copyProb,
								parentCutSet)));
						//Create a holding queue to hold results for expansion after the callable queue shrinks.
						LinkedBlockingQueue<Tuple2<double[], HashSet<Integer>>> holdingQueue = new
								LinkedBlockingQueue<Tuple2<double[], HashSet<Integer>>>();
						//Take jobs off the front and either enqueue derivative problems or push their result to a
						//queue for expansion once the callable queue size is manageable again.
						while (!callableQueue.isEmpty() || !holdingQueue.isEmpty()) {
							//If the callableQueue size drops below 10, and there are solutions to be expanded, we
							//want to expand a stored solution.
							if (callableQueue.size() < 10 && !holdingQueue.isEmpty()) {
								expandSolution(holdingQueue.poll(), callableQueue, lpArray.get(i), completedCutSets,
										pool);
							}
							//Pull a result off the queue:
							Tuple2<double[], HashSet<Integer>> resultTuple = callableQueue.poll().get();
							//Record the result:
							int total = 0;
							HashSet<Integer> activeEntitiesSet = new HashSet<Integer>();
							for (int index = 0; index < resultTuple.a.length; index++) {
								if (objectiveSet.contains(index) && resultTuple.a[index] > 0.0) {
									total++;
								}
								if (resultTuple.a[index] > 0.0) {
									activeEntitiesSet.add(index);
								}
							}
							Tuple2<HashSet<Integer>, HashSet<Integer>> representation
									= new Tuple2<HashSet<Integer>, HashSet<Integer>>(activeEntitiesSet, resultTuple.b);
							if (!resultMap.containsKey(total)) {
								ArrayList<Tuple2<HashSet<Integer>, HashSet<Integer>>> al
										= new ArrayList<Tuple2<HashSet<Integer>, HashSet<Integer>>>();
								al.add(representation);
								resultMap.put(total, al);
							} else {
								resultMap.get(total).add(representation);
							}
							//Push the result onto the holding queue for subsequent cutset generation:
							holdingQueue.put(resultTuple);
						}
						TextOutputHandler.writeResultConditionComparison(this, resultMap, i);
					}
				}
			} else {
				log.error("INVALID PROBLEM TYPE " + this.mode);
				System.exit(-1);
			}
		} catch (Exception e) {
			log.error("ERROR", e);
			e.printStackTrace();
		} finally {
			if (pool != null) {
				pool.shutdownNow();
			}
			GLPK.glp_delete_prob(mainLP);
			mainLP=null;
		}
	}

	private void expandSolution(Tuple2<double[], HashSet<Integer>> resultTuple,
								LinkedBlockingQueue<Future<Tuple2<double[], HashSet<Integer>>>> callableQueue,
								glp_prob lpSolveInstance, HashSet<HashSet<Integer>> completedCutSets, ExecutorService
										pool)
			{
		//Unpack result from Tuple
		double[] resultArray = resultTuple.a;
		HashSet<Integer> cutSet = resultTuple.b;
		if (resultArray.length == 1) {
			//No result - terminate expansion (no solution found, etc).
			return;
		} else {
			//Have result - generate derivative cuts.
			for (Integer index : this.objectiveSet) {
				if (resultArray[index] > 0.0) {
					HashSet<Integer> newCutSet = new HashSet<Integer>(cutSet);
					//Cut this value from the result.
					newCutSet.add(index);
					//Check to see if we've already processed this cutset:
					if (completedCutSets.contains(newCutSet)) {
						continue;
					} else {
						//Instantiate a new LP and enqueue it:
						glp_prob newProb = GLPK.glp_create_prob();
						GLPK.glp_copy_prob(newProb, lpSolveInstance, GLPK.GLP_ON);
						IntegerProgrammingTools.cut(newProb, newCutSet);
						callableQueue.add(pool.submit(new LpSolveCutSetCallable(newProb, newCutSet)));
						completedCutSets.add(newCutSet);
					}
				}
			}
		}

	}

	public File getOutputDirectory() {
		return outputDirectory;
	}

	public HashMap<Object, Integer> getSetterMap() {
		return this.setterMap;
	}

	public boolean getSaveLowConf() {
		return saveLowConf;
	}

	public boolean getSaveMidConf() {
		return saveMidConf;
	}

	public String getDescriptor() {
		return descriptor;
	}
}
