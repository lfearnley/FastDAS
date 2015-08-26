package ModelCuration;

import FastDAS.KeyParam;
import LogicalSystem.Tuple2;
import LogicalSystem.Util;
import org.apache.log4j.Logger;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.*;
import org.biopax.paxtools.model.level3.Process;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.LinkedBlockingQueue;

/**
 * Created with IntelliJ IDEA.
 * User: uqlfearn
 * Date: 26/08/13
 * Time: 12:56 PM
 * To change this template use File | Settings | File Templates.
 */
public class BucketFixer {

	private static org.apache.log4j.Logger log = Logger.getLogger(BucketFixer.class);

	private static HashMap<PhysicalEntity, Integer> countMap;


	public static boolean needsExpansion(PhysicalEntity pe) {
		boolean expand = false;
		//Two possibilities - either a Complex containing EntitySet components somewhere:
		if (pe instanceof Complex && !((Complex) pe).getComponent().isEmpty()) {
			//Step through each component
			for (PhysicalEntity component : ((Complex) pe).getComponent()) {
				//If that component needs expansion, return true:
				if (needsExpansion(component)) {
					expand = true;
					break;
				}
			}
			//Or it's a straight EntitySet, in which case we check members for disjoint reaction involvement.
		} else if (pe.getComment().contains(KeyParam.REACTOME_ENTITYSET_FLAG)) {
			//Run through each member physicalEntity:
			for (PhysicalEntity member : pe.getMemberPhysicalEntity()) {
				//If the member participates in reactions independently of its entity set then we need to flag for
				// expansion
				if (!member.getParticipantOf().isEmpty()) {
					expand = true;
					break;
				}
			}
		}
		return expand;
	}

	private static Tuple2<HashSet<PhysicalEntity>, HashSet<Interaction>> extractSubnetwork(PhysicalEntity pe) {
		HashSet<Interaction> iset = new HashSet<Interaction>();
		HashSet<PhysicalEntity> peset = new HashSet<PhysicalEntity>();
		Tuple2<HashSet<PhysicalEntity>, HashSet<Interaction>> subNetTuple
				= new Tuple2<HashSet<PhysicalEntity>, HashSet<Interaction>>(peset, iset);
		try {
			LinkedBlockingQueue<Entity> entityQueue = new LinkedBlockingQueue<Entity>();
			entityQueue.put(pe);
			while (!entityQueue.isEmpty()) {
				Entity current = entityQueue.poll();
				if (current instanceof PhysicalEntity &&
						(current.getComment().contains(KeyParam.REACTOME_COMPLEX_ENTITYSET_CONTAINING_FLAG) ||
								current.getComment().contains(KeyParam.REACTOME_ENTITYSET_FLAG))) {
					//TODO - code to test whether the entityset is derivative of PE or not....
					peset.add((PhysicalEntity) current);
					for (Interaction i : current.getParticipantOf()) {
						if (!iset.contains(i)) {
							entityQueue.add(i);
						}
					}
				} else if (current instanceof Control) {
					iset.add((Interaction) current);
				} else if (current instanceof Conversion) {
					iset.add((Interaction) current);
					for (PhysicalEntity leftPE : ((Conversion) current).getLeft()) {
						if (!peset.contains(leftPE)) {
							entityQueue.add(leftPE);
						}
					}
					for (PhysicalEntity rightPE : ((Conversion) current).getRight()) {
						if (!peset.contains(rightPE)) {
							entityQueue.add(rightPE);
						}
					}
				}
			}
			return subNetTuple;
		} catch (InterruptedException e) {
			e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
		}
		return subNetTuple;
	}

	public static class PriorityComparator<PhysicalEntity> implements Comparator<PhysicalEntity> {
		private HashMap<PhysicalEntity, Integer> peMap;

		public PriorityComparator(HashMap<PhysicalEntity, Integer> peMap) {
			this.peMap = peMap;
		}

		public int compare(PhysicalEntity n1, PhysicalEntity n2) {
			return countMap.get(n1).compareTo(countMap.get(n2));
		}
	}

	private static boolean isEquivalent(PhysicalEntity pe1, PhysicalEntity pe2) {
		if (pe1.equals(pe2)) {
			return true;
		} else if (!pe1.getCellularLocation().equals(pe2.getCellularLocation()) && pe1.getDisplayName().equals(pe2
				.getDisplayName())) {
			return true;
		} else if (pe1.getComment().contains(KeyParam.REACTOME_ENTITYSET_FLAG) && pe2.getComment().contains(KeyParam
				.REACTOME_ENTITYSET_FLAG)) {
			for (PhysicalEntity memberPE : pe1.getMemberPhysicalEntity()) {
				boolean map = false;
				for (PhysicalEntity memberPE2 : pe2.getMemberPhysicalEntity()) {
					if (isEquivalent(memberPE, memberPE2)) {
						map = true;
						break;
					}
				}
				if (!map) {
					return false;
				}
			}
			return true;
		} else if (pe1 instanceof Protein && pe2 instanceof Protein) {
			//Capture protein modification by checking entity references.
			if (((Protein) pe1).getEntityReference() != null && ((Protein) pe2).getEntityReference() != null) {
				return ((Protein) pe1).getEntityReference().equals(((Protein) pe2).getEntityReference());
			} else {
				return false;
			}
		} else if (pe1 instanceof Complex && pe2 instanceof Complex) {
			Complex complex1 = (Complex) pe1;
			Complex complex2 = (Complex) pe2;
			//Complex assembly case
			if (complex1.getComponent().contains(complex2) || complex2.getComponent().contains(complex1)) {
				return true;
				//} else if (Util.hasComplexOverlap(complex1, complex2) || Util.hasComplexOverlap(complex2,
				// complex1)) {
				//    return true;
			} else if (complexMaps(complex1, complex2)) {
				return true;
			}
		} else if (pe1 instanceof Complex || pe2 instanceof Complex) {
			Complex complex = (pe1 instanceof Complex) ? (Complex) pe1 : (Complex) pe2;
			PhysicalEntity entity = !(pe1 instanceof Complex) ? pe1 : pe2;
			if (complex.getComponent().contains(entity)) {
				return true;
			} else if (Util.hasEquivalentInComplex(complex, entity)) {
				return true;
			}
		} else {
			log.warn("Issue in equivalency testing pe1 : " + pe1.getDisplayName() + " " + pe1.getRDFId() + " to " +
					pe2.getDisplayName() + " " + pe2.getRDFId());
			return false;
		}
		return false;
	}

	private static boolean complexMaps(Complex complex1, Complex complex2) {
		boolean matchedforward = true;
		boolean matchedreverse = true;
		for (PhysicalEntity pe : complex1.getComponent()) {
			boolean matched = false;
			if (complex2.getComponent().contains(pe)) {
				continue;
			}
			for (PhysicalEntity pe2 : complex2.getComponent()) {
				if (pe instanceof Protein && pe2 instanceof Protein && !pe.getComment().contains(KeyParam
						.REACTOME_ENTITYSET_FLAG) && !pe2.getComment().contains(KeyParam.REACTOME_ENTITYSET_FLAG)) {
					if (((Protein) pe).getEntityReference().equals(((Protein) pe2).getEntityReference())) {
						matched = true;
						break;
					}
				} else if (pe instanceof Protein && pe2 instanceof Protein) {
					return pe.equals(pe2);
				}
			}
			if (!matched) {
				matchedforward = false;
				break;
			}
		}
		if (matchedforward) {
			return true;
		}
		for (PhysicalEntity pe : complex2.getComponent()) {
			boolean matched = false;
			if (complex1.getComponent().contains(pe)) {
				continue;
			}
			for (PhysicalEntity pe2 : complex1.getComponent()) {
				if (pe instanceof Protein && pe2 instanceof Protein && !pe.getComment().contains(KeyParam
						.REACTOME_ENTITYSET_FLAG) && !pe2.getComment().contains(KeyParam.REACTOME_ENTITYSET_FLAG)) {
					if (((Protein) pe).getEntityReference().equals(((Protein) pe2).getEntityReference())) {
						matched = true;
						break;
					}
				} else if (pe instanceof Protein && pe2 instanceof Protein) {
					return pe.equals(pe2);
				}
			}
			if (!matched) {
				matchedreverse = false;
				break;
			}
		}
		return matchedreverse;
	}

	public static void debucketModel(Model model, File outputDirectory) {
		try {
			BufferedWriter logWriter = new BufferedWriter(new FileWriter(new File(outputDirectory.getAbsolutePath() +
					System.getProperty("file.separator") + "debucketerLog.txt")));
			long expansionCount = 0;
			HashSet<PhysicalEntity> peset = new HashSet<PhysicalEntity>(model.getObjects(PhysicalEntity.class));
			int lim = 50;
			int count = 0;
			for (PhysicalEntity pe : peset) {
				if (!pe.getComment().contains(KeyParam.REACTOME_ENTITYSET_FLAG)) {
					//Test to see whether this is in an EntitySet.
					HashMap<PhysicalEntity, PhysicalEntity> peMap = new HashMap<PhysicalEntity, PhysicalEntity>();
					if (!pe.getParticipantOf().isEmpty() && !pe.getMemberPhysicalEntityOf().isEmpty()) {
						int escount = 0;
						for (PhysicalEntity es : pe.getMemberPhysicalEntityOf()) {
							escount += es.getParticipantOf().size();
						}
						if (escount == 0) {
							continue;
						}
						logWriter.write("****************************************");
						logWriter.newLine();
						count++;
						logWriter.write("RUNNING " + pe.getDisplayName() + "(" + pe.getRDFId() + ") in " + pe
								.getParticipantOf().size() + " reactions " + pe.getMemberPhysicalEntityOf().size() +
								"" +
								" entitySets in " + escount);
						HashSet<Interaction> done = new HashSet<Interaction>();
						for (PhysicalEntity entitySet : pe.getMemberPhysicalEntityOf()) {
							logWriter.write("PROCESSING : " + entitySet.getDisplayName() + " (" + entitySet.getRDFId()
									+ ")");
							logWriter.newLine();
							peMap.put(entitySet, pe);
							LinkedBlockingQueue<Interaction> convQueue = new LinkedBlockingQueue<Interaction>();
							for (Interaction i : entitySet.getParticipantOf()) {
								if (!done.contains(i)) {
									done.add(i);
									convQueue.add(i);
									break;
								}
							}
							exit:
							while (!convQueue.isEmpty()) {
								Interaction currentInteraction = convQueue.poll();
								if (currentInteraction instanceof Conversion) {
									Conversion currentConversion = (Conversion) currentInteraction;
									logWriter.write("CURRENTLY ON " + currentConversion.getDisplayName() +
											currentConversion.getRDFId());
									logWriter.newLine();
									for (PhysicalEntity lpe : currentConversion.getLeft()) {
										logWriter.write(lpe.getDisplayName() + " (" + lpe.getRDFId().replaceAll(model
												.getXmlBase(), "") + ") + ");
									}
									logWriter.write("----------->");
									for (PhysicalEntity rpe : currentConversion.getRight()) {
										logWriter.write(rpe.getDisplayName() + " (" + rpe.getRDFId().replaceAll(model
												.getXmlBase(), "") + ") + ");
									}
									logWriter.newLine();

									boolean ok = false;
									out:
									for (PhysicalEntity leftHandSideEntity : currentConversion.getLeft()) {
										if (peMap.containsKey(leftHandSideEntity)) {
											logWriter.write("FOUND " + leftHandSideEntity.getDisplayName());
											logWriter.newLine();
											//Util.replaceInConversion(newConv, leftHandSideEntity, peMap.get
											// (leftHandSideEntity));
											for (PhysicalEntity rightHandSideEntity : currentConversion.getRight()) {
												if (isEquivalent(leftHandSideEntity, rightHandSideEntity)) {
													logWriter.write("WE MAPPED " + leftHandSideEntity.getDisplayName()
															+ " TO " + rightHandSideEntity.getDisplayName());
													logWriter.newLine();
													logWriter.write("WE REPLACED LEFT " + leftHandSideEntity
															.getDisplayName() + " WITH " + peMap.get
															(leftHandSideEntity).getDisplayName());
													logWriter.newLine();
													if (rightHandSideEntity.isEquivalent(leftHandSideEntity)) {
														Conversion newConv = Util.duplicateConversion
																(currentConversion, model, expansionCount++, model
																		.getXmlBase());
														Util.replaceInConversion(newConv, leftHandSideEntity, peMap
																.get(leftHandSideEntity));
														Util.replaceInConversion(newConv, rightHandSideEntity, peMap
																.get(leftHandSideEntity));
														ok = true;
														logWriter.write("WAS EQUAL");
														logWriter.newLine();
													}
													if (!leftHandSideEntity.getCellularLocation().equals
															(rightHandSideEntity.getCellularLocation()) &&
															leftHandSideEntity.getDisplayName().equals
																	(rightHandSideEntity.getDisplayName())) {
														//We want to duplicate the replacement that we had
														PhysicalEntity replacement = Util.duplicateWithTransport
																(model, peMap.get(leftHandSideEntity),
																		rightHandSideEntity, expansionCount);
														logWriter.write("TRANSPORT: AND RIGHT " + replacement
																.getDisplayName() + " (" + replacement.getRDFId() +
																")" +
																" in " + rightHandSideEntity.getParticipantOf().size
																());
														logWriter.newLine();
														peMap.put(rightHandSideEntity, replacement);
														Conversion newConv = Util.duplicateConversion
																(currentConversion, model, expansionCount++, model
																		.getXmlBase());
														Util.replaceInConversion(newConv, leftHandSideEntity, peMap
																.get(leftHandSideEntity));
														Util.replaceInConversion(newConv, rightHandSideEntity, peMap
																.get(rightHandSideEntity));
														ok = true;
														for (Interaction nextInteraction : rightHandSideEntity
																.getParticipantOf()) {
															if (!done.contains(nextInteraction)) {
																done.add(nextInteraction);
																convQueue.add(nextInteraction);
															}
														}
														expansionCount++;
													} else if (rightHandSideEntity instanceof Complex &&
															rightHandSideEntity.getMemberPhysicalEntity().isEmpty()) {
														logWriter.write("INTO HERE");
														for (PhysicalEntity component : ((Complex)
																rightHandSideEntity).getComponent()) {
															if (peMap.containsKey(component)) {
																Complex newComplex = Util.getEquivalentComplex(model,
																		(Complex) rightHandSideEntity, component,
																		peMap.get(component), expansionCount);
																ok = true;
																expansionCount++;
																logWriter.write("AND RIGHT " + rightHandSideEntity
																		.getDisplayName() + " WITH " + newComplex
																		.getDisplayName());
																logWriter.newLine();
																peMap.put(rightHandSideEntity, newComplex);
																Conversion newConv = Util.duplicateConversion
																		(currentConversion, model, expansionCount++,
																				model.getXmlBase());
																Util.replaceInConversion(newConv, leftHandSideEntity,
																		peMap.get(leftHandSideEntity));
																Util.replaceInConversion(newConv, rightHandSideEntity,
																		peMap.get(rightHandSideEntity));
																for (Interaction nextInteraction : rightHandSideEntity
																		.getParticipantOf()) {
																	if (!done.contains(nextInteraction)) {
																		done.add(nextInteraction);
																		convQueue.add(nextInteraction);
																	}
																}
																break out;
															}
														}
													} else {
														for (PhysicalEntity replacement : getReplacementEntities
																(rightHandSideEntity)) {
															logWriter.write("TRYING : " + replacement.getDisplayName
																	());
															logWriter.newLine();
															if (isEquivalent(peMap.get(leftHandSideEntity),
																	replacement)) {
																logWriter.write("AND RIGHT " + rightHandSideEntity
																		.getDisplayName() + " WITH " + replacement
																		.getDisplayName());
																logWriter.newLine();
																ok = true;
																peMap.put(rightHandSideEntity, replacement);
																Conversion newConv = Util.duplicateConversion
																		(currentConversion, model, expansionCount++,
																				model.getXmlBase());
																Util.replaceInConversion(newConv, leftHandSideEntity,
																		peMap.get(leftHandSideEntity));
																Util.replaceInConversion(newConv, rightHandSideEntity,
																		peMap.get(rightHandSideEntity));
																for (Interaction nextInteraction : rightHandSideEntity
																		.getParticipantOf()) {
																	if (!done.contains(nextInteraction)) {
																		done.add(nextInteraction);
																		convQueue.add(nextInteraction);
																	}
																}
																break out;
															}
														}
													}
												}
												if (!ok) {
													logWriter.write("FAILURE");
													logWriter.newLine();
													logWriter.write(leftHandSideEntity.getDisplayName());
													logWriter.newLine();
													logWriter.write(peMap.get(leftHandSideEntity).getDisplayName());
													logWriter.newLine();
													if (peMap.get(leftHandSideEntity) instanceof Complex) {
														for (PhysicalEntity component : ((Complex) peMap.get
																(leftHandSideEntity)).getComponent()) {
															logWriter.write("\t" + component.getDisplayName() + " (" +
																	component.getRDFId() + ")");
															logWriter.newLine();
														}
													}
													logWriter.newLine();
													break exit;
												}
											}
										}
									}
									if (!ok) {
										out2:
										for (PhysicalEntity rightHandSideEntity : currentConversion.getRight()) {
											if (peMap.containsKey(rightHandSideEntity)) {
												logWriter.write("RIGHT LOOP FOUND " + rightHandSideEntity
														.getDisplayName());
												logWriter.newLine();
												//Util.replaceInConversion(newConv, leftHandSideEntity, peMap.get
												// (leftHandSideEntity));
												for (PhysicalEntity leftHandSideEntity : currentConversion.getLeft()) {
													if (isEquivalent(rightHandSideEntity, leftHandSideEntity)) {
														logWriter.write("WE MAPPED " + rightHandSideEntity
																.getDisplayName() + " TO " + rightHandSideEntity
																.getDisplayName());
														logWriter.newLine();
														logWriter.write("WE REPLACED RIGHT " + rightHandSideEntity
																.getDisplayName() + " WITH " + peMap.get
																(rightHandSideEntity).getDisplayName());
														logWriter.newLine();
														if (leftHandSideEntity.isEquivalent(rightHandSideEntity)) {
															ok = true;
															Conversion newConv = Util.duplicateConversion
																	(currentConversion, model, expansionCount++, model
																			.getXmlBase());
															Util.replaceInConversion(newConv, leftHandSideEntity,
																	peMap.get(leftHandSideEntity));
															Util.replaceInConversion(newConv, rightHandSideEntity,
																	peMap.get(rightHandSideEntity));
															logWriter.write("WAS EQUAL");
															logWriter.newLine();
														}
														if (!rightHandSideEntity.getCellularLocation().equals
																(leftHandSideEntity.getCellularLocation()) &&
																rightHandSideEntity.getDisplayName().equals
																		(leftHandSideEntity.getDisplayName())) {
															//We want to duplicate the replacement that we had
															PhysicalEntity replacement = Util.duplicateWithTransport
																	(model, peMap.get(rightHandSideEntity),
																			leftHandSideEntity, expansionCount);
															logWriter.write("TRANSPORT: AND LEFT " + replacement
																	.getDisplayName() + " (" + replacement.getRDFId()
																	+ ") in " + leftHandSideEntity.getParticipantOf()
																	.size());
															logWriter.newLine();
															peMap.put(leftHandSideEntity, replacement);
															Conversion newConv = Util.duplicateConversion
																	(currentConversion, model, expansionCount++, model
																			.getXmlBase());
															Util.replaceInConversion(newConv, leftHandSideEntity,
																	peMap.get(leftHandSideEntity));
															Util.replaceInConversion(newConv, rightHandSideEntity,
																	peMap.get(rightHandSideEntity));
															ok = true;
															for (Interaction nextInteraction : leftHandSideEntity
																	.getParticipantOf()) {
																if (!done.contains(nextInteraction)) {
																	done.add(nextInteraction);
																	convQueue.add(nextInteraction);
																}
															}
															expansionCount++;
														} else if (leftHandSideEntity instanceof Complex &&
																leftHandSideEntity.getMemberPhysicalEntity().isEmpty
																		()) {
															logWriter.write("INTO HERE");
															for (PhysicalEntity component : ((Complex)
																	leftHandSideEntity).getComponent()) {
																if (peMap.containsKey(component)) {
																	Complex newComplex = Util.getEquivalentComplex
																			(model, (Complex) leftHandSideEntity,
																					component, peMap.get(component),
																					expansionCount);
																	ok = true;
																	expansionCount++;
																	logWriter.write("AND LEFT " + leftHandSideEntity
																			.getDisplayName() + " WITH " + newComplex
																			.getDisplayName());
																	logWriter.newLine();
																	peMap.put(leftHandSideEntity, newComplex);
																	Conversion newConv = Util.duplicateConversion
																			(currentConversion, model,
																					expansionCount++, model.getXmlBase
																							());
																	Util.replaceInConversion(newConv,
																			leftHandSideEntity, peMap.get
																					(leftHandSideEntity));
																	Util.replaceInConversion(newConv,
																			rightHandSideEntity, peMap.get
																					(rightHandSideEntity));
																	for (Interaction nextInteraction :
																			leftHandSideEntity.getParticipantOf()) {
																		if (!done.contains(nextInteraction)) {
																			done.add(nextInteraction);
																			convQueue.add(nextInteraction);
																		}
																	}
																	break out2;
																}
															}
														} else {
															for (PhysicalEntity replacement : getReplacementEntities
																	(leftHandSideEntity)) {
																logWriter.write("TRYING : " + replacement
																		.getDisplayName());
																logWriter.newLine();
																if (isEquivalent(peMap.get(rightHandSideEntity),
																		replacement)) {
																	logWriter.write("AND LEFT " + leftHandSideEntity
																			.getDisplayName() + " WITH " + replacement
																			.getDisplayName());
																	logWriter.newLine();
																	ok = true;
																	peMap.put(leftHandSideEntity, replacement);
																	Conversion newConv = Util.duplicateConversion
																			(currentConversion, model,
																					expansionCount++, model.getXmlBase
																							());
																	Util.replaceInConversion(newConv,
																			leftHandSideEntity, peMap.get
																					(leftHandSideEntity));
																	Util.replaceInConversion(newConv,
																			rightHandSideEntity, peMap.get
																					(rightHandSideEntity));
																	for (Interaction nextInteraction :
																			leftHandSideEntity.getParticipantOf()) {
																		if (!done.contains(nextInteraction)) {
																			done.add(nextInteraction);
																			convQueue.add(nextInteraction);
																		}
																	}
																	break out2;
																}
															}
														}
													}
													if (!ok) {
														logWriter.write("FAILURE");
														logWriter.newLine();
														logWriter.write(rightHandSideEntity.getDisplayName());
														logWriter.newLine();
														logWriter.write(peMap.get(rightHandSideEntity).getDisplayName
																());
														logWriter.newLine();
														if (peMap.get(rightHandSideEntity) instanceof Complex) {
															for (PhysicalEntity component : ((Complex) peMap.get
																	(rightHandSideEntity)).getComponent()) {
																logWriter.write("\t" + component.getDisplayName() + "" +
																		" " +
																		"(" + component.getRDFId() + ")");
																logWriter.newLine();
															}
														}
														logWriter.newLine();
														break exit;
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		} catch (Exception e) {
			log.error("Exception in debucketer.", e);
			System.exit(-1);
		}
	}

	private static HashSet<PhysicalEntity> getReplacementEntities(PhysicalEntity pe) {
		HashSet<PhysicalEntity> replacements = new HashSet<PhysicalEntity>();
		//Run through each member physicalEntity:
		for (PhysicalEntity member : pe.getMemberPhysicalEntity()) {
			//If the member participates in reactions independently of its entity set then we need to flag for
			// expansion
			if (!member.getParticipantOf().isEmpty()) {
				replacements.add(member);
			}
		}
		return replacements;
	}

	public static void writeEntitySetNetwork(Model model) {
		try {
			BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(new File("out.psi")));
			BufferedWriter bufferedattrWriter = new BufferedWriter(new FileWriter(new File("out.noa")));
			HashSet<Interaction> interactionHashSet = new HashSet<Interaction>();
			HashSet<PhysicalEntity> peset = new HashSet<PhysicalEntity>();
			HashSet<PhysicalEntity> written = new HashSet<PhysicalEntity>();
			for (PhysicalEntity pe : model.getObjects(PhysicalEntity.class)) {
				if (pe.getComment().contains(KeyParam.REACTOME_ENTITYSET_FLAG) || pe.getComment().contains(KeyParam
						.REACTOME_COMPLEX_ENTITYSET_CONTAINING_FLAG)) {
					if (needsExpansion(pe)) {
						peset.add(pe);
						interactionHashSet.addAll(pe.getParticipantOf());
					}
				}
			}
			for (Interaction i : interactionHashSet) {
				if (i instanceof Conversion) {
					Conversion conv = ((Conversion) i);
					String convName = "";
					if (conv.getDisplayName() == null) {
						convName = conv.getRDFId().replaceAll(model.getXmlBase(), "");
						log.warn("WARN : " + conv.getRDFId().replaceAll(model.getXmlBase(), "") + " has no " +
								"displayname");
					} else {
						convName = conv.getDisplayName().replaceAll(" ", "_");
					}
					for (PhysicalEntity pe : conv.getLeft()) {
						if (peset.contains(pe)) {
							written.add(pe);
							String peName = "";
							if (pe.getDisplayName() == null) {
								peName = pe.getRDFId().replaceAll(model.getXmlBase(), "");
								log.warn("WARN : " + pe.getRDFId().replaceAll(model.getXmlBase(), "") + " has no " +
										"displayname");
							} else {
								peName = pe.getDisplayName().replaceAll(" ", "_");
							}
							bufferedWriter.write(peName + "\tpp\t" + convName);
							bufferedWriter.newLine();
						}
					}
					for (PhysicalEntity pe : conv.getRight()) {
						if (peset.contains(pe)) {
							written.add(pe);
							String peName = "";
							if (pe.getDisplayName() == null) {
								peName = pe.getRDFId().replaceAll(model.getXmlBase(), "");
								log.warn("WARN : " + pe.getRDFId().replaceAll(model.getXmlBase(), "") + " has no " +
										"displayname");
							} else {
								peName = pe.getDisplayName().replaceAll(" ", "_");
							}
							bufferedWriter.write(convName + "\tpp\t" + peName);
							bufferedWriter.newLine();
						}
					}
				} else if (i instanceof Control) {
					Control control = (Control) i;
					HashSet<String> processSet = new HashSet<String>();
					for (Process proc : control.getControlled()) {
						if (proc instanceof Conversion) {
							String procName = "";
							if (proc.getDisplayName() == null) {
								procName = "CONTROL_" + proc.getRDFId().replaceAll(model.getXmlBase(), "");
								log.warn("WARN : " + proc.getRDFId().replaceAll(model.getXmlBase(), "") + " has no " +
										"displayname");
							} else {
								procName = proc.getDisplayName().replaceAll(" ", "_");
							}
							processSet.add(procName);
						}
					}
					HashSet<String> controllersSet = new HashSet<String>();
					for (Controller ctrlr : control.getController()) {
						if (peset.contains(ctrlr)) {
							String peName = "";
							if (ctrlr.getDisplayName() == null) {
								peName = ctrlr.getRDFId().replaceAll(model.getXmlBase(), "");
								log.warn("WARN : " + ctrlr.getRDFId().replaceAll(model.getXmlBase(), "") + " has no " +
										"displayname");
							} else {
								peName = ctrlr.getDisplayName().replaceAll(" ", "_");
							}
							written.add((PhysicalEntity) ctrlr);
							controllersSet.add(peName);
						}
					}
					for (String controller : controllersSet) {
						for (String process : processSet) {
							bufferedWriter.write(controller + "\tctrl\t" + process);
							bufferedWriter.newLine();
						}
					}
				}
			}
			bufferedWriter.flush();
			bufferedWriter.close();
			for (PhysicalEntity pe : written) {
				boolean isinput = true;
				for (Interaction i : pe.getParticipantOf()) {
					if (i instanceof Conversion) {
						Conversion conv = (Conversion) i;
						if (conv.getRight().contains(pe)) {
							isinput = false;
						}
					}
				}
				String peName = "";
				if (pe.getDisplayName() == null) {
					peName = pe.getRDFId().replaceAll(model.getXmlBase(), "");
					log.warn("WARN : " + pe.getRDFId().replaceAll(model.getXmlBase(), "") + " has no displayname");
				} else {
					peName = pe.getDisplayName().replaceAll(" ", "_");
				}
				bufferedattrWriter.write(peName + "\t" + isinput);
				bufferedattrWriter.newLine();

			}
			bufferedattrWriter.flush();
			bufferedattrWriter.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}


}
