package IO;

import FastDAS.KeyParam;
import LogicalSystem.Util;
import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.*;
import org.biopax.paxtools.model.level3.Process;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.LinkedBlockingQueue;

/**
 * Created with IntelliJ IDEA.
 * User: uqlfearn
 * Date: 15/08/13
 * Time: 3:56 PM
 * To change this template use File | Settings | File Templates.
 */
public class BioPAXInputHandler {

	private static org.apache.log4j.Logger log = Logger.getLogger(BioPAXInputHandler.class);

	/**
	 * Gets a BioPAX model from a file passed as an argument, then processes it for use with FastDAS.
	 *
	 * @param biopaxFile
	 * @param pathwayFile
	 * @param speciesFile
	 * @return
	 */
	public static Model getModel(File biopaxFile, File pathwayFile, File speciesFile) {
		Model model = readModelFromFile(biopaxFile);
		//Uncomment below line to dump a plain-text representation of the pathway hierarchy.
		//TextOutputHandler.writePathwayHierarchy(model, new File("PathwayHierarchy.txt"));
		flattenPathways(model);
		//Strip non-biochemical reactions from the model.
		reduceModel(model);
		//TODO: Check
		HashSet<Control> ctrlremoves = new HashSet<Control>();
		for (Conversion conv : model.getObjects(Conversion.class)) {
			for (Control ctrl : conv.getControlledOf()) {
				HashSet<Controller> removes = new HashSet<Controller>();
				for (Controller ctrlr : ctrl.getController()) {
					if (ctrlr instanceof PhysicalEntity && conv.getLeft().contains(ctrlr)) {
						removes.add(ctrlr);
					}
				}
				if (!removes.isEmpty()) {
					for (Controller ctrlr : removes) {
						ctrl.removeController(ctrlr);
					}
					if (ctrl.getController().isEmpty()) {
						ctrlremoves.add(ctrl);
					}
				}
			}
		}
		for (Control ctrl : ctrlremoves) {
			model.remove(ctrl);
		}
		//Strip out pathways for deletion.
		if (pathwayFile==null) {
			log.warn("WARN: No pathway file specified.");
		} else {
			removePathways(pathwayFile, model);
		}
		//Strip out species for deletion.
		if (speciesFile==null) {
			log.warn("WARN: No species file specified.");
		} else {
			removeSpecies(speciesFile, model);
			//Populate experimental conditions specified in the species file.
			populateExperimentalConditions(speciesFile, model);
		}
		//Perform a sanity check - look for nonsense in the reactions and interactions in the file.
		performSanityCheck(model);
		return model;
	}

	/**
	 * Method for stripping out any interaction that isn't a Conversion or a Control. Necessary to limit to purely
	 * biochemical reaction (i.e. stripping out DNA synthesis, RNA synthesis, etc).
	 *
	 * @param model A paxtools Model object.
	 */
	private static void reduceModel(Model model) {
		HashSet<Interaction> iset = new HashSet<Interaction>(model.getObjects(Interaction.class));
		for (Interaction currentInteraction : iset) {
			if (!(currentInteraction instanceof Conversion || currentInteraction instanceof Control)) {
				Util.removeInteraction(currentInteraction, model);
			}
		}
	}

	/**
	 * Brief sanity check using some basic heuristics for model validation. This checks to see whether there are any
	 * reactions without substrates or products, or any control events without an appropriate controller. Reactions
	 * and control events meeting these criteria are removed from the model and logged. The method also checks for
	 * missing cellular locations, and will log these but take no action.
	 *
	 * @param model A paxtools Model object.
	 */
	private static void performSanityCheck(Model model) {
		HashSet<Interaction> interactionHashSet = new HashSet<Interaction>(model.getObjects(Interaction.class));
		for (Interaction currentInteraction : interactionHashSet) {
			if (currentInteraction instanceof Conversion) {
				Conversion currentConversion = (Conversion) currentInteraction;
				if (currentConversion.getLeft().isEmpty() || currentConversion.getRight().isEmpty()) {
					log.warn("WARN : " + currentConversion.getRDFId().replaceAll(model.getXmlBase(), "") + " has a " +
							"missing product or substrate.");
					model.remove(currentConversion);
				}
			} else if (currentInteraction instanceof Control) {
				Control currentControl = (Control) currentInteraction;
				if (currentControl.getControlled().isEmpty() || currentControl.getController().isEmpty()) {
					model.remove(currentControl);
				} else if (currentControl.getControlType() == null) {
					log.warn("WARN : " + currentControl.getRDFId().replaceAll(model.getXmlBase(), "") + " has no " +
							"controlType set");
					model.remove(currentControl);
				}
			} else if (currentInteraction instanceof TemplateReaction) {
				TemplateReaction currentTemplateReaction = (TemplateReaction) currentInteraction;
				if (currentTemplateReaction.getProduct().isEmpty()) {
					model.remove(currentTemplateReaction);
				}
			}
		}
		//Check and log missing cellular locations - database errors, etc.
		for (PhysicalEntity pe : model.getObjects(PhysicalEntity.class)) {
			if (pe.getCellularLocation() == null) {
				log.warn("WARN : " + pe.getRDFId().replaceAll(model.getXmlBase(), "") + " has no cellular location");
			}
		}
	}

	/**
	 * Flatten pathways, removing any ordering from the reactions. Pathway objects contain interactions ordered by
	 * pathwaysteps or defined via subpathways. We want to merge these for our purposes (we do not assume that any
	 * ordering is imposed (as ordering makes little sense given the amount of crosstalk that can occur within the
	 * system).
	 *
	 * @param m A paxtools model object.
	 */
	public static void flattenPathways(Model m) {
		HashSet<Pathway> pathset = new HashSet<Pathway>(m.getObjects(Pathway.class));
		//Iterate over the pathways in the model
		for (Pathway p : pathset) {
			LinkedBlockingQueue<PathwayStep> processQueue = new LinkedBlockingQueue<PathwayStep>();
			processQueue.addAll(p.getPathwayOrder());
			while (!processQueue.isEmpty()) {
				PathwayStep currentProcess = processQueue.poll();
				for (Process childProcess : (currentProcess.getStepProcess())) {
					if (childProcess instanceof Interaction) {
						p.addPathwayComponent(childProcess);
					}
				}
				p.removePathwayOrder(currentProcess);
				m.remove(currentProcess);
			}
		}
	}

	/**
	 * Removes all species marked for deletion in a FastDAS species file from the model.  Designed to log deletion
	 * failures rather than throwing an
	 * exception specifically to
	 *
	 * @param speciesFile A FastDAS format species file.
	 * @param m           A paxtools model object.
	 */
	public static void removeSpecies(File speciesFile, Model m) {
		try {
			BufferedReader speciesFileReader = new BufferedReader(new FileReader(speciesFile));
			String currentLine = null;
			while ((currentLine = speciesFileReader.readLine()) != null) {
				String[] lineArray = currentLine.split("\t");
				if (lineArray[KeyParam.SPECIES_DELETION_FLAG].trim().equals("REMOVE")) {
					if (m.getByID(m.getXmlBase() + lineArray[1]) == null) {
						continue;
					} else {
						PhysicalEntity pe = (PhysicalEntity) m.getByID(m.getXmlBase() + lineArray[KeyParam
								.SPECIES_RDF_IDX]);
						HashSet<Interaction> iset = new HashSet<Interaction>(pe.getParticipantOf());
						for (Interaction i : iset) {
							if (i instanceof Conversion) {
								((Conversion) i).removeLeft(pe);
								((Conversion) i).removeRight(pe);
							} else if (i instanceof Control) {
								((Control) i).removeController(pe);
							} else {
								log.warn("Attempting to remove " + lineArray[KeyParam.SPECIES_RDF_IDX] + " (" +
										lineArray[KeyParam.SPECIES_NAME_IDX] + "), participating in " + i.getRDFId());
								continue;
							}
						}
					}
				}
			}
		} catch (Exception e) {
			log.error("Error loading species for deletion", e);
			System.exit(-1);
		}
	}

	/**
	 * This method takes a FastDAS species file and populates experimental conditions by adding comments to species.
	 * If activity is in the CONTROL condition, this is tagged with "CONTROL ", otherwise "EXP_<number>" followed by
	 * the condition number. Specific flags for activity are ACTIVE/INACTIVE/FREE.
	 *
	 * @param speciesFile A FastDAS format species file.
	 * @param m           A paxtools model object.
	 */
	public static void populateExperimentalConditions(File speciesFile, Model m) {
		try {
			//Open species file.
			BufferedReader speciesFileReader = new BufferedReader(new FileReader(speciesFile));
			String currentLine;
			int conditionCount = 1;
			while ((currentLine = speciesFileReader.readLine()) != null) {
				String[] lineArray = currentLine.split("\t");
				//Check to see if there are any specified conditions.
				if (lineArray.length == KeyParam.SPECIES_CONDITION_IDX_START) {
					continue;
				}
				//Check to see if the RDF is still in the model...
				if (m.getByID(m.getXmlBase() + lineArray[KeyParam.SPECIES_RDF_IDX]) == null) {
					//If not, then continue; log this for the user.
					//log.warn("WARN : " + lineArray[KeyParam.SPECIES_RDF_IDX] + " is no longer in model and cannot " +
							//"have experimental conditions set.");
					continue;
				} else {
					//Found species, set its activity.
					PhysicalEntity pe = (PhysicalEntity) m.getByID(m.getXmlBase() + lineArray[KeyParam
							.SPECIES_RDF_IDX]);
					//Set control activity.
					pe.addComment("CONTROL " + lineArray[KeyParam.SPECIES_CONDITION_IDX_START]);
					//For each experimental condition set by the user...
					int idx = KeyParam.SPECIES_CONDITION_IDX_START + 1;
					while (idx < lineArray.length) {
						//Add that experimental condition to the species.
						String comment = "EXP_" + (idx - KeyParam.SPECIES_CONDITION_IDX_START) + lineArray[idx].trim();
						pe.addComment(comment);
						idx++;
					}
					//TEST!!!!
					/*HashSet<Interaction> iset = new HashSet<Interaction>();
					iset = new HashSet<Interaction>(pe.getParticipantOf());
					for (Interaction i : iset) {
						if (i instanceof Conversion) {
							Conversion conv = (Conversion)i;
							if (conv.getRight().contains(pe)) {
								conv.removeRight(pe);
							}
						}
					}
					*/
					//Update condition counter if necessary.
					if (conditionCount < (idx - KeyParam.SPECIES_CONDITION_IDX_START))  {
						System.out.println("INCREMENT : " + lineArray[1] + " WAS " + conditionCount);
						conditionCount = (idx - KeyParam.SPECIES_CONDITION_IDX_START);
						System.out.println(" NOW " + conditionCount);
					}
				}
			}
			//Update conditioncount object within the model.
			//TODO: Check to see if there's a model property that would be better to use than this kludge.
			PhysicalEntity cc = m.addNew(PhysicalEntity.class, "ConditionCount");
			cc.addComment("" + conditionCount);
			System.out.println("SETTING " + conditionCount);
		} catch (Exception e) {
			log.error("Error loading species for deletion", e);
			System.exit(-1);
		}
	}

	/**
	 * This method removes a pathway, its sub-pathways, and any reactions annotated as being unique members of this
	 * pathway. We do not delete a reaction if it is a member of more than one pathway - in this case, we remove the
	 * annotation for the deleted pathway from the reaction.
	 *
	 * @param pathwayFile A FastDAS formatted pathway removal file.
	 * @param model       A paxtools model object.
	 */
	public static void removePathways(File pathwayFile, Model model) {
		HashSet<Process> pathsToRemove = getPathwaysForRemoval(pathwayFile, model);
		HashSet<Interaction> processesToRemove = new HashSet<Interaction>();
		HashSet<Interaction> interactionSet = new HashSet<Interaction>(model.getObjects(Interaction.class));
		for (Interaction currentInteraction : interactionSet) {
			Set<Pathway> pathwaySet = currentInteraction.getPathwayComponentOf();
			if (pathwaySet.isEmpty()) {
				continue;
			}
			boolean remove = true;
			for (Pathway p : pathwaySet) {
				if (!pathsToRemove.contains(p)) {
					remove = false;
				}
			}
			if (remove) {
				processesToRemove.add(currentInteraction);
			}
		}
		HashSet<Entity> participantSet = new HashSet<Entity>(model.getObjects(Entity.class));
		HashSet<Entity> participantsToRemove = new HashSet<Entity>();
		for (Entity participant : participantSet) {
			Set<Interaction> participationSet = participant.getParticipantOf();
			if (participationSet.isEmpty()) {
				continue;
			}
			boolean remove = true;
			for (Interaction interaction : participationSet) {
				if (!processesToRemove.contains(interaction)) {
					remove = false;
					break;
				}
			}
			if (remove) {
				participantsToRemove.add(participant);
			}
		}
		for (Interaction i : processesToRemove) {
			model.remove(i);
		}
		for (Process p : pathsToRemove) {
			model.remove(p);
		}
	}

	/**
	 * Helper method - takes a list of pathways, and generates a list of its subpathways and processes in order to
	 * remove them in BioPAXInputHandler.removePathways().
	 *
	 * @param pathwayFile A FastDAS formatted pathway removal file.
	 * @param m           A paxtools model object.
	 * @return
	 */
	private static HashSet<Process> getPathwaysForRemoval(File pathwayFile, Model m) {
		HashSet<Process> pathwaySet = new HashSet<Process>();
		try {
			BufferedReader pathwayReader = new BufferedReader(new FileReader(pathwayFile));
			String current;
			LinkedBlockingQueue<Pathway> pathwayQueue = new LinkedBlockingQueue<Pathway>();
			//Line by line processing of pathway file
			while ((current = pathwayReader.readLine()) != null) {
				Pathway pathway = getPathwayByName(current.trim(), m);
				if (pathway == null) {
					continue;
				}
				pathwayQueue.add(pathway);
				pathwaySet.add(pathway);
			}
			while (!pathwayQueue.isEmpty()) {
				Pathway currentpath = pathwayQueue.poll();
				for (Process p : currentpath.getPathwayComponent()) {
					if (p instanceof Pathway && !pathwaySet.contains(p)) {
						pathwayQueue.add((Pathway) p);
						pathwaySet.add(p);
					}
				}
			}
			return pathwaySet;
		} catch (Exception e) {
			//Failure case - return an empty pathway set, warn in log.
			log.warn("Pathway set failed to load - skipping pathway reduction step");
			return new HashSet<Process>();
		}
	}

	/**
	 * Helper method that gets pathway objects by their display name. Does not support duplicate names.
	 * TODO: Sort out duplicate pathway names in this method.
	 *
	 * @param pathwayName a String containing the display name of a Pathway object
	 * @param m           A paxtools model object
	 * @return the Pathway with the given display name
	 */
	private static Pathway getPathwayByName(String pathwayName, Model m) {
		//Iterate through all pathways in the model.
		for (Pathway p : m.getObjects(Pathway.class)) {
			//If the pathway names include the target name, return it.
			if (p.getName().contains(pathwayName)) {
				return p;
			}
		}
		//Otherwise, warn that the pathway is not contained in the database and return a null result.
		log.warn("Could not find pathway named " + pathwayName);
		return null;
	}

	/**
	 * Reads a BioPAX Level 3 file into memory, returning a paxtools Model object.
	 *
	 * @param biopaxFile a java.io.File instance containing BioPAX Level 3 data (e.g. new File("homo sapiens.owl"))
	 * @return a paxtools Model object specified in the provided file
	 */
	public static Model readModelFromFile(File biopaxFile) {
		FileInputStream fis = null;
		Model model = null;
		try {
			SimpleIOHandler ioHandler = new SimpleIOHandler();
			fis = new FileInputStream(biopaxFile);
			model = ioHandler.convertFromOWL(fis);
		} catch (Exception e) {
			if (biopaxFile == null) {
				log.error("Read of BioPAX model failed - null file, check your file selection/file contents.");
				log.error("Shutting down.");
				System.exit(-1);
			} else {
				String filename = biopaxFile.getName();
				String errorstring = "Read of BioPAX model failed - check " + filename;
				log.error(errorstring, e);
				log.error("Shutting down.");
				System.exit(-1);
			}
		} finally {
			if (fis != null) {
				IOUtils.closeQuietly(fis);
			}
			return model;
		}
	}
}
