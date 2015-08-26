package IO;

import FastDAS.KeyParam;
import LogicalSystem.AdjList;
import LogicalSystem.ProblemInstance;
import LogicalSystem.Util;
import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;
import org.biopax.paxtools.model.level3.*;
import org.biopax.paxtools.model.level3.Process;
import sun.awt.image.ImageWatched;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.security.Key;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;

/**
 * Created with IntelliJ IDEA.
 * User: uqlfearn
 * Date: 16/08/13
 * Time: 1:04 PM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkOutputHandler {

	private static org.apache.log4j.Logger log = Logger.getLogger(NetworkOutputHandler.class);

	private static ArrayList<HashSet<Integer>> upstreamReactionTraverser(ProblemInstance instance, double[] results) {
		HashSet<Integer> specifiedSpecies = new HashSet<Integer>();
		HashMap<Object, Integer> setterMap = instance.getSetterMap();
		HashMap<Integer, Object> integerObjectHashMap = instance.getIntegerObjectHashMap();
		HashMap<Object, Integer> objectIntegerHashMap = instance.getObjectIntegerHashMap();
		for (Object o : setterMap.keySet()) {
			Integer settercode = setterMap.get(o);
			if (settercode == KeyParam.SETTER_USER_SPEC || settercode == KeyParam.SETTER_ESET_FILL || settercode ==
					KeyParam.SETTER_COMPLEX_FILL) {
				specifiedSpecies.add(objectIntegerHashMap.get(o));
			}
		}
		HashSet<Integer> upstreamActiveReactions = new HashSet<Integer>();
		HashSet<Integer> downstreamActiveReactions = new HashSet<Integer>();
		AdjList al = instance.getNetworkAdjList();
		HashSet<Integer> seen = new HashSet<Integer>();
		LinkedBlockingQueue<Integer> traversalQueue = new LinkedBlockingQueue<Integer>();
		traversalQueue.addAll(specifiedSpecies);
		seen.addAll(specifiedSpecies);
		while (!traversalQueue.isEmpty()) {
			Integer currentInteger = traversalQueue.poll();
			HashSet<Integer> parents = al.getParents(Math.abs(currentInteger));
			for (Integer i : parents) {
				if (!seen.contains(Math.abs(i))) {
					seen.add(Math.abs(i));
					traversalQueue.add(Math.abs(i));
					if ((integerObjectHashMap.get(Math.abs(i))) instanceof Interaction
							&& results[i] != 0.0) {
						upstreamActiveReactions.add(Math.abs(i));
					}
				}
			}
		}
		traversalQueue.clear();
		seen.clear();
		seen.addAll(specifiedSpecies);
		traversalQueue.addAll(specifiedSpecies);
		while (!traversalQueue.isEmpty()) {
			Integer currentInteger = traversalQueue.poll();
			HashSet<Integer> children = al.getChildren(Math.abs(currentInteger));
			for (Integer i : children) {
				if (!seen.contains(Math.abs(i))) {
					seen.add(Math.abs(i));
					traversalQueue.add(Math.abs(i));
					if ((integerObjectHashMap.get(Math.abs(i))) instanceof Interaction
							&& results[i] != 0.0) {
						downstreamActiveReactions.add(Math.abs(i));
					}
				}
			}
		}
		System.out.println(upstreamActiveReactions.size() + "  " + downstreamActiveReactions.size());
		HashSet<Integer> internalActiveReactions = new HashSet<Integer>(upstreamActiveReactions);
		internalActiveReactions.retainAll(downstreamActiveReactions);
		upstreamActiveReactions.removeAll(internalActiveReactions);
		downstreamActiveReactions.removeAll(internalActiveReactions);
		System.out.println("THERE WERE " + specifiedSpecies.size() + " specified species joined by " +
				internalActiveReactions.size() + " with " + upstreamActiveReactions.size() +" upstream reactions and " +
				downstreamActiveReactions.size() +" downstream reactions");
		ArrayList<HashSet<Integer>> returnList = new ArrayList<HashSet<Integer>>();
		returnList.add(internalActiveReactions);
		returnList.add(upstreamActiveReactions);
		returnList.add(downstreamActiveReactions);
		return returnList;
	}

	/**
	 * Method for writing a differential graph and associated files (node attributes, affected proteins, affected
	 * pathways, and downstream pathways/proteins from the solver system.
	 *
	 * @param results              ArrayList of Future double arrays from the GPLK solver queue
	 * @param instance             ProblemInstance object describing the problem being formulated
	 * @param allResultsInOneGraph flag - if true, all results are written to a single graph rather than pairwise
	 *                             comparisons
	 */
	public static void writeDifferentialGraph(ArrayList<Future<double[]>> results, ProblemInstance instance, boolean
			allResultsInOneGraph) {
		try {
			//Create output folder(s) as required.
			File outputFolder = instance.getOutputDirectory();
			if (!instance.getDescriptor().equals("")) {
				outputFolder = new File(instance.getOutputDirectory().getPath() + KeyParam.PATH_SEPARATOR + instance
						.getDescriptor());
			}
			boolean success = outputFolder.mkdirs();
			File highQualOutputFolder = new File(outputFolder.getPath() + KeyParam.PATH_SEPARATOR + "High_Conf");
			File midQualOutputFolder = new File(outputFolder.getPath() + KeyParam.PATH_SEPARATOR + "Mid_Conf");
			File lowQualOutputFolder = new File(outputFolder.getPath() + KeyParam.PATH_SEPARATOR + "Low_Conf");
			File unionGraphOutputFolder = new File(outputFolder.getPath() + KeyParam.PATH_SEPARATOR + "Union_Graph");
			highQualOutputFolder.mkdirs();
			midQualOutputFolder.mkdirs();
			lowQualOutputFolder.mkdirs();
			unionGraphOutputFolder.mkdirs();
			//Generate high quality prediction lists to output.
			//Ordered list (indexes corresponding to experiment indexes) of interactions to write out.
			ArrayList<HashSet<Integer>> highQualityPredictionSetList = new ArrayList<HashSet<Integer>>();
			ArrayList<HashSet<Integer>> mediumQualityPredictionSetList = new ArrayList<HashSet<Integer>>();
			ArrayList<HashSet<Integer>> lowQualityPredictionSetList = new ArrayList<HashSet<Integer>>();
			//Get the result for the 'control' experiment.
			double[] controlResult = results.get(0).get();
			//For each experimental simulation (i.e., from condition '1' onwards):
			for (int i = 1; i < results.size(); i++) {
				//Difference set initiation.
				HashSet<Integer> differenceSet = new HashSet<Integer>();
				double[] currentResult = results.get(i).get();
				//For each interaction in the problem
				System.out.println("REACTION INDICES\n----------------");
				for (Integer reactionIndex : instance.getReactionSet()) {
					System.out.println(reactionIndex);
					//Test for difference.
					if (controlResult[reactionIndex] != currentResult[reactionIndex]) {
						System.out.println("GRABBING : " + reactionIndex);
						differenceSet.add(reactionIndex);
					}
				}
				//Add the difference set to the set for writing out.
				highQualityPredictionSetList.add(differenceSet);
			}
			//TODO: Write out high quality here.
			writeDiffGraph(highQualityPredictionSetList, instance, allResultsInOneGraph, highQualOutputFolder);
			//Generate mid-confidence prediction lists for output.
			if (instance.getSaveMidConf()) {
				for (int i = 0; i < highQualityPredictionSetList.size(); i++) {
					HashSet<Integer> medQualityPredictions = generatePredictionOutputSets
							(highQualityPredictionSetList.get(i), instance, true);
					medQualityPredictions.removeAll(highQualityPredictionSetList.get(i));
					mediumQualityPredictionSetList.add(medQualityPredictions);
				}
				//TODO: Write out mid quality here.
				writeDiffGraph(mediumQualityPredictionSetList, instance, allResultsInOneGraph, midQualOutputFolder);
			}
			//Generate low-confidence prediction lists for output.
			if (instance.getSaveLowConf()) {
				for (int i = 0; i < highQualityPredictionSetList.size(); i++) {
					HashSet<Integer> lowQualityPredictions = generatePredictionOutputSets
							(highQualityPredictionSetList.get(i), instance, false);
					lowQualityPredictions.removeAll(highQualityPredictionSetList.get(i));
					lowQualityPredictions.removeAll(mediumQualityPredictionSetList.get(i));
					lowQualityPredictionSetList.add(lowQualityPredictions);
				}
				//TODO: Write out low quality here.
				writeDiffGraph(lowQualityPredictionSetList, instance, allResultsInOneGraph, lowQualOutputFolder);
			}
			if (instance.getSaveLowConf() || instance.getSaveMidConf()) {
				//TODO: Write out union network and predictions.
				writeQualityMergedDiffGraph(highQualityPredictionSetList, mediumQualityPredictionSetList,
						lowQualityPredictionSetList, instance, allResultsInOneGraph,
						unionGraphOutputFolder);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static void writeQualityMergedDiffGraph(ArrayList<HashSet<Integer>> highSetList,
													ArrayList<HashSet<Integer>> medSetList,
													ArrayList<HashSet<Integer>> lowSetList,
													ProblemInstance instance, boolean allResultsInOneGraph,
													File outputFolder) {
		String fileName = new SimpleDateFormat("yyyyMMddhhmm").format(new Date());
		//Create BufferedWriters - GML, Node Attribute, affected Pathway, and affected Protein
		BufferedWriter bufferedGMLWriter = null;
		BufferedWriter bufferedNodeAttributeWriter = null;
		BufferedWriter bufferedPathwayWriter = null;
		BufferedWriter bufferedProteinWriter = null;
		try {
			ArrayList<HashSet<Integer>> outputsetlist = new ArrayList<HashSet<Integer>>();
			for (int i = 0; i < highSetList.size(); i++) {
				HashSet<Integer> currentSet = new HashSet<Integer>();
				currentSet.addAll(highSetList.get(i));
				currentSet.addAll(medSetList.get(i));
				currentSet.addAll(lowSetList.get(i));
				outputsetlist.add(currentSet);
			}
			//Instantiate the Node Attribute writer.
			bufferedNodeAttributeWriter = new BufferedWriter(new FileWriter(new File(outputFolder.getPath()
					+ System.getProperty("file.separator") + "NodeAttributes" + fileName + ".txt")));
			//Get a set of all interactions to be written (union of all differences).
			HashSet<Integer> unionInteractionsSet = new HashSet<Integer>();
			for (HashSet<Integer> hs : outputsetlist) {
				unionInteractionsSet.addAll(hs);
			}
			//If we're writing all results in one graph (NOT RECOMMENDED):
			if (allResultsInOneGraph) {
				//START WRITER INSTANTIATION
				bufferedGMLWriter = new BufferedWriter(new FileWriter(new File(outputFolder.getPath()
						+ System.getProperty("file.separator") + "Network_" + "combined_" + fileName + ".gml")));
				//END WRITER INSTANTIATION
				writeGMLPreamble(bufferedGMLWriter, instance, fileName, "combined");
				writeGMLFile(bufferedGMLWriter, instance, unionInteractionsSet, fileName, "combined");
				writePathwayFile(bufferedPathwayWriter, instance, unionInteractionsSet, fileName, "combined", outputFolder);
				writeProteinList(bufferedProteinWriter, instance, unionInteractionsSet, fileName, "combined", outputFolder);
			//If we're writing results in separate graphs:
			} else {
				int instanceCount = 1;
				for (HashSet<Integer> hashSet : outputsetlist) {
					//START WRITER INSTANTIATION
					bufferedGMLWriter = new BufferedWriter(new FileWriter(new File(outputFolder.getPath()
							+ System.getProperty("file.separator") + "Network_" + "Control-" + instanceCount + "_" +
							fileName + ".gml")));
					//END WRITER INSTANTIATION
					writeGMLPreamble(bufferedGMLWriter, instance, fileName, "Control-" + instanceCount);
					writeGMLFile(bufferedGMLWriter, instance, hashSet, fileName, "Control-" + instanceCount);
					writePathwayFile(bufferedPathwayWriter, instance, hashSet, fileName, "Control-" + instanceCount,
							outputFolder);
					writeProteinList(bufferedProteinWriter, instance, hashSet, fileName, "Control-" + instanceCount,
							outputFolder);
					instanceCount++;
				}
			}
			//Write the node attribute file.
			writeQualMergedNodeAttributeFile(bufferedNodeAttributeWriter, instance, unionInteractionsSet,
					highSetList, medSetList, lowSetList);
			//Pokemon handling - quite a bit can go wrong here, so we're catching it all.
		} catch (Exception e) {
			log.error("ERROR : Differential writer error", e);
			//Close all output streams, gracefully close.
		} finally {
			IOUtils.closeQuietly(bufferedGMLWriter);
			IOUtils.closeQuietly(bufferedNodeAttributeWriter);
			IOUtils.closeQuietly(bufferedPathwayWriter);
			IOUtils.closeQuietly(bufferedProteinWriter);
		}
	}

	private static void writeDiffGraph(ArrayList<HashSet<Integer>> outputsetlist, ProblemInstance instance, boolean
			allResultsInOneGraph, File outputFolder) {
		//TODO: Filename is currently a date/time string. This probably needs an update.
		String fileName = new SimpleDateFormat("yyyyMMddhhmm").format(new Date());
		//Create BufferedWriters - GML, Node Attribute, affected Pathway, and affected Protein
		BufferedWriter bufferedGMLWriter = null;
		BufferedWriter bufferedNodeAttributeWriter = null;
		BufferedWriter bufferedPathwayWriter = null;
		BufferedWriter bufferedProteinWriter = null;
		try {
			//Instantiate the Node Attribute writer.
			bufferedNodeAttributeWriter = new BufferedWriter(new FileWriter(new File(outputFolder.getPath()
					+ System.getProperty("file.separator") + "NodeAttributes" + fileName + ".txt")));
			//Get a set of all interactions to be written (union of all differences).
			HashSet<Integer> unionInteractionsSet = new HashSet<Integer>();
			for (HashSet<Integer> hs : outputsetlist) {
				unionInteractionsSet.addAll(hs);
			}
			//If we're writing all results in one graph (NOT RECOMMENDED):
			if (allResultsInOneGraph) {
				//START WRITER INSTANTIATION
				bufferedGMLWriter = new BufferedWriter(new FileWriter(new File(outputFolder.getPath()
						+ System.getProperty("file.separator") + "Network_" + "combined_" + fileName + ".gml")));
				//END WRITER INSTANTIATION
				writeGMLPreamble(bufferedGMLWriter, instance, fileName, "combined");
				writeGMLFile(bufferedGMLWriter, instance, unionInteractionsSet, fileName, "combined");
				writePathwayFile(bufferedPathwayWriter, instance, unionInteractionsSet, fileName, "combined",
						outputFolder);
				writeProteinList(bufferedProteinWriter, instance, unionInteractionsSet, fileName, "combined",
						outputFolder);
				//If we're writing results in separate graphs:
			} else {
				int instanceCount = 1;
				for (HashSet<Integer> hashSet : outputsetlist) {
					//START WRITER INSTANTIATION
					bufferedGMLWriter = new BufferedWriter(new FileWriter(new File(outputFolder.getPath()
							+ System.getProperty("file.separator") + "Network_" + "Control-" + instanceCount + "_" +
							fileName + ".gml")));
					//END WRITER INSTANTIATION
					writeGMLPreamble(bufferedGMLWriter, instance, fileName, "Control-" + instanceCount);
					writeGMLFile(bufferedGMLWriter, instance, hashSet, fileName, "Control-" + instanceCount);
					writePathwayFile(bufferedPathwayWriter, instance, hashSet, fileName, "Control-" + instanceCount,
							outputFolder);
					writeProteinList(bufferedProteinWriter, instance, hashSet, fileName, "Control-" + instanceCount,
							outputFolder);
					instanceCount++;
				}
			}
			//Write the node attribute file.
			writeNodeAttributeFile(bufferedNodeAttributeWriter, instance, unionInteractionsSet, outputFolder);
			//Pokemon handling - quite a bit can go wrong here, so we're catching it all.
		} catch (Exception e) {
			log.error("ERROR : Differential writer error", e);
			//Close all output streams, gracefully close.
		} finally {
			IOUtils.closeQuietly(bufferedGMLWriter);
			IOUtils.closeQuietly(bufferedNodeAttributeWriter);
			IOUtils.closeQuietly(bufferedPathwayWriter);
			IOUtils.closeQuietly(bufferedProteinWriter);
		}
	}

	private static void writeQualMergedNodeAttributeFile(BufferedWriter bufferedNodeAttributeWriter, ProblemInstance
			instance, HashSet<Integer> interactionSet, ArrayList<HashSet<Integer>> highQualList,
														 ArrayList<HashSet<Integer>> medQualList,
														 ArrayList<HashSet<Integer>> lowQualList)
			throws IOException {
		HashSet<Integer> linesToWrite = new HashSet<Integer>();
		AdjList adjList = instance.getNetworkAdjList();
		for (Integer i : interactionSet) {
			linesToWrite.add(Math.abs(i));
			linesToWrite.addAll(adjList.getNeighbours(Math.abs(i)));
		}
		HashMap<Integer, Object> integerObjectHashMap = instance.getIntegerObjectHashMap();
		for (Integer currentInteger : linesToWrite) {
			StringBuilder stringBuilder = new StringBuilder();
			Entity currentEntity = (Entity) integerObjectHashMap.get(currentInteger);
			if (currentEntity == null) {
				continue;
			}
			stringBuilder.append(currentInteger);
			stringBuilder.append("\t");
			String name = getNodeName(currentEntity);
			stringBuilder.append(name);
			stringBuilder.append("\t");
			stringBuilder.append(currentEntity.getRDFId().replaceAll(instance.getModel().getXmlBase(), ""));
			stringBuilder.append("\t");
			//Add cellular location if PhysicalEntity, otherwise skip (Interactions don't have locations):
			if (currentEntity instanceof Interaction) {
				stringBuilder.append("\t");
			} else {
				CellularLocationVocabulary cellularLocationVocabulary = ((PhysicalEntity) currentEntity)
						.getCellularLocation();
				if (cellularLocationVocabulary != null) {
					String locationString = cellularLocationVocabulary.getTerm().iterator().next();
					if (locationString == null) {
						stringBuilder.append(cellularLocationVocabulary.toString());
					} else {
						stringBuilder.append(locationString);
					}
				} else {
					//Logger warns in getNodeName method if there's no cellular location associated.
					stringBuilder.append("Not Recorded");
				}
				stringBuilder.append("\t");
			}
			//Add URL link to DB to node:
			if (currentEntity instanceof Interaction) {
				if (currentEntity.getXref() != null) {
					HashSet<Xref> xrefs = new HashSet<Xref>(currentEntity.getXref());
					for (Xref xr : xrefs) {
						if (xr.getDb()!=null && xr.getDb().equals("Reactome")) {
							stringBuilder.append(KeyParam.REACT_URL);
							stringBuilder.append(xr.getId());
						}
					}
				}
			} else {
				stringBuilder.append("\t");
			}
			if (!(currentEntity instanceof Interaction) && instance.getSetterMap().containsKey(currentEntity)) {
				int settercode = instance.getSetterMap().get(currentEntity);
				if (settercode == KeyParam.SETTER_USER_SPEC) {
					stringBuilder.append("USER");
				} else if (settercode == KeyParam.SETTER_COMPLEX_FILL) {
					stringBuilder.append("COMPLEX FILL");
				} else if (settercode == KeyParam.SETTER_ESET_FILL) {
					stringBuilder.append("ESET FILL");
				}
			} else {
				stringBuilder.append("\t");
			}
			for (int i = 0; i < highQualList.size(); i++) {
				if (highQualList.get(i).contains(currentInteger)) {
					stringBuilder.append("Exp_" + (i + 1) + " HIGH\t");
				} else if (medQualList.get(i).contains(currentInteger)) {
					stringBuilder.append("Exp_" + (i + 1) + " MEDIUM\t");
				} else if (lowQualList.get(i).contains(currentInteger)) {
					stringBuilder.append("Exp_" + (i + 1) + " LOW\t");
				} else {
					stringBuilder.append("\t");
				}
			}
			stringBuilder.append(System.getProperty("line.separator"));
			bufferedNodeAttributeWriter.write(stringBuilder.toString());
		}
		bufferedNodeAttributeWriter.flush();
		bufferedNodeAttributeWriter.close();
	}

	private static void writeProteinList(BufferedWriter bufferedProteinWriter, ProblemInstance instance,
										 HashSet<Integer> interactionIndexSet, String fileName, String modifier, File
										 outputFolder) {
		try {
			bufferedProteinWriter = new BufferedWriter(new FileWriter(new File(outputFolder.getPath() +
					System.getProperty("file.separator") + "Proteins" + fileName + "_" + modifier + ".txt")));
			HashMap<Integer, Object> integerObjectHashMap = instance.getIntegerObjectHashMap();
			HashSet<Interaction> interactionSet = new HashSet<Interaction>();
			for (Integer interactionIndex : interactionIndexSet) {
				interactionSet.add((Interaction) integerObjectHashMap.get(interactionIndex));
			}
			HashSet<Entity> doneEntitySet = new HashSet<Entity>();
			for (Interaction interaction : interactionSet) {
				LinkedBlockingQueue<Entity> entityLinkedBlockingQueue = new LinkedBlockingQueue<Entity>();
				entityLinkedBlockingQueue.addAll(interaction.getParticipant());
				while (!entityLinkedBlockingQueue.isEmpty()) {
					Entity entity = entityLinkedBlockingQueue.poll();
					if (!doneEntitySet.contains(entity)) {
						if (entity instanceof Protein) {
							entityLinkedBlockingQueue.addAll(((Protein) entity).getMemberPhysicalEntity());
							Protein protein = (Protein) entity;
							bufferedProteinWriter.write(getNodeName(protein));
							bufferedProteinWriter.write("\t");
							bufferedProteinWriter.write(Util.getUniProt(protein));
							bufferedProteinWriter.write("\t");
							bufferedProteinWriter.newLine();
						} else if (entity instanceof PhysicalEntity && !((PhysicalEntity) entity)
								.getMemberPhysicalEntity().isEmpty()) {
							entityLinkedBlockingQueue.addAll(((PhysicalEntity) entity).getMemberPhysicalEntity());
						} else if (entity instanceof Complex) {
							entityLinkedBlockingQueue.addAll(((Complex) entity).getComponent());
						}
						doneEntitySet.add(entity);
					}
				}
			}
			bufferedProteinWriter.flush();
		} catch (Exception e) {
			log.error("ERROR: Writing Pathway file " + fileName + modifier + " threw an error", e);
		}
	}

	private static void writePathwayFile(BufferedWriter bufferedPathwayWriter, ProblemInstance instance,
										 HashSet<Integer> interactionIndexSet, String fileName, String modifier, File
										 outputFolder) {
		try {
			bufferedPathwayWriter = new BufferedWriter(new FileWriter(new File(outputFolder.getPath() +
					KeyParam.PATH_SEPARATOR + "Pathway" + fileName + "_" + modifier + ".txt")));
			HashMap<Pathway, Integer> pathwayCountMap = new HashMap<Pathway, Integer>();
			HashMap<Integer, Object> integerObjectHashMap = instance.getIntegerObjectHashMap();
			HashSet<Interaction> interactionSet = new HashSet<Interaction>();
			for (Integer interactionIndex : interactionIndexSet) {
				interactionSet.add((Interaction) integerObjectHashMap.get(Math.abs(interactionIndex)));
			}
			for (Interaction interaction : interactionSet) {
				for (Pathway path : interaction.getPathwayComponentOf()) {
					if (!pathwayCountMap.containsKey(path)) {
						pathwayCountMap.put(path, 1);
					} else {
						pathwayCountMap.put(path, pathwayCountMap.get(path) + 1);
					}
				}
			}
			for (Pathway pathway : pathwayCountMap.keySet()) {
				bufferedPathwayWriter.write(pathway.getDisplayName() + "\t" + pathwayCountMap.get(pathway));
				bufferedPathwayWriter.newLine();
			}
			bufferedPathwayWriter.flush();
		} catch (Exception e) {
			log.error("ERROR: Writing Pathway file " + fileName + modifier + " threw an error", e);
		}
	}

	private static void writeNodeAttributeFile(BufferedWriter bufferedNodeAttributeWriter, ProblemInstance instance,
											   HashSet<Integer> interactionSet, File outputFolder) throws IOException {
		HashSet<Integer> linesToWrite = new HashSet<Integer>();
		AdjList adjList = instance.getNetworkAdjList();
		for (Integer i : interactionSet) {
			linesToWrite.add(Math.abs(i));
			linesToWrite.addAll(adjList.getNeighbours(Math.abs(i)));
		}
		bufferedNodeAttributeWriter.write("ID\tDisplay_Name\tRDF_ID\tTYPE\tCELLULAR_LOCATION\tURL\tSETTER\tINHIBITOR" +
				KeyParam
						.NEWLINE);
		HashMap<Integer, Object> integerObjectHashMap = instance.getIntegerObjectHashMap();
		for (Integer currentInteger : linesToWrite) {
			StringBuilder stringBuilder = new StringBuilder();
			Entity currentEntity = (Entity) integerObjectHashMap.get(Math.abs(currentInteger));
			if (currentEntity == null) {
				continue;
			}
			stringBuilder.append(currentInteger);
			stringBuilder.append("\t");
			String name = getNodeName(currentEntity);
			stringBuilder.append(name);
			stringBuilder.append("\t");
			stringBuilder.append(currentEntity.getRDFId().replaceAll(instance.getModel().getXmlBase(), ""));
			stringBuilder.append("\t");
			if (currentEntity instanceof Interaction) {
				stringBuilder.append("INTERACTION\t");
			} else if (currentEntity instanceof Protein) {
				stringBuilder.append("PROTEIN\t");
			} else if (currentEntity instanceof Complex) {
				stringBuilder.append("COMPLEX\t");
			} else if (currentEntity instanceof Control) {
				stringBuilder.append("CONTROL\t");
			} else if (currentEntity instanceof SmallMolecule) {
				stringBuilder.append("SMALL_MOLECULE\t");
			} else if (currentEntity instanceof PhysicalEntity) {
				stringBuilder.append("PHYSICAL_ENTITY\t");
			} else {
				stringBuilder.append("OUTPUT_ERROR_UNCAUGHT_TYPE\t");
			}
			//Add cellular location if PhysicalEntity, otherwise skip (Interactions don't have locations):
			if (currentEntity instanceof Interaction) {
				stringBuilder.append("\t");
			} else {
				CellularLocationVocabulary cellularLocationVocabulary = ((PhysicalEntity) currentEntity)
						.getCellularLocation();
				if (cellularLocationVocabulary != null) {
					String locationString = cellularLocationVocabulary.getTerm().iterator().next();
					if (locationString == null) {
						stringBuilder.append(cellularLocationVocabulary.toString());
					} else {
						stringBuilder.append(locationString);
					}
				} else {
					//Logger warns in getNodeName method if there's no cellular location associated.
					stringBuilder.append("Not Recorded");
				}
				stringBuilder.append("\t");
			}
			//Add URL link to DB to node:
			if (currentEntity instanceof Interaction) {
				if (currentEntity.getXref() != null) {
					HashSet<Xref> xrefs = new HashSet<Xref>(currentEntity.getXref());
					for (Xref xr : xrefs) {
						if (xr.getDb() !=null && xr.getDb().equals("Reactome")) {
							stringBuilder.append(KeyParam.REACT_URL);
							stringBuilder.append(xr.getId());
						}
					}
				}
			} else {
				stringBuilder.append("\t");
			}
			if (!(currentEntity instanceof Interaction) && instance.getSetterMap().containsKey(currentEntity)) {
				int settercode = instance.getSetterMap().get(currentEntity);
				if (settercode == KeyParam.SETTER_USER_SPEC) {
					stringBuilder.append("USER");
				} else if (settercode == KeyParam.SETTER_COMPLEX_FILL) {
					stringBuilder.append("COMPLEX FILL");
				} else if (settercode == KeyParam.SETTER_ESET_FILL) {
					stringBuilder.append("ESET FILL");
				}
			}
			stringBuilder.append("\t");
			if (currentInteger < 0) {
				stringBuilder.append("INHIBITOR");
			}
			stringBuilder.append(System.getProperty("line.separator"));
			bufferedNodeAttributeWriter.write(stringBuilder.toString());
		}
		bufferedNodeAttributeWriter.flush();
		bufferedNodeAttributeWriter.close();
	}

	private static String getNodeName(Entity currentEntity) {
		try {
			String name = "";
			//Name behaviour - default to displayName for the entity. This should always be set. If not, then go to
			// standard
			//name, then any other recorded name. If the entity has no names associated, then spit out the RDF ID -
			// this is
			//considered a database error and will trigger a log warning.
			name = currentEntity.getDisplayName();
			if (name == null) {
				name = currentEntity.getStandardName();
				if (name == null) {
					if (currentEntity.getName() != null && !currentEntity.getName().isEmpty()) {
						name = currentEntity.getName().iterator().next();
						if (name == null) {
							name = currentEntity.getRDFId();
							log.warn(currentEntity.getRDFId() + " has no associated name.");
						}
					} else {
						name = currentEntity.getRDFId();
						log.warn(currentEntity.getRDFId() + " has no associated name.");
					}
				}
			}
			//Cellular location - add to end of name if the entity is a PhysicalEntity.
			String cellLocation = "";
			if (currentEntity instanceof PhysicalEntity) {
				CellularLocationVocabulary cellularLocationVocabulary = ((PhysicalEntity) currentEntity)
						.getCellularLocation();
				//Error condition - no cellular location being recorded.
				if (cellularLocationVocabulary == null) {
					//Log as warning.
					log.warn(currentEntity.getRDFId() + " has no associated location");
				} else {
					String locationString = cellularLocationVocabulary.getTerm().iterator().next();
					if (locationString == null) {
						locationString = cellularLocationVocabulary.toString();
					}
					cellLocation = " " + locationString;
				}
			}
			name = name + cellLocation;
			return name;
		} catch (Exception e) {
			log.error("EXCEPTION IN GET NAME : " + currentEntity.getRDFId(), e);
			return currentEntity.getRDFId();
		}
	}

	private static void writeGMLFile(BufferedWriter bufferedGMLWriter, ProblemInstance instance, HashSet<Integer>
			interactionSet,
									 String fileName, String nameModifierString) throws IOException {
		AdjList adjList = instance.getNetworkAdjList();
		HashSet<Integer> written = new HashSet<Integer>();
		StringBuilder sb = new StringBuilder();
		for (Integer i : interactionSet) {
			if (!written.contains(Math.abs(i))) {
				written.add(Math.abs(i));
				bufferedGMLWriter.write("\tnode [");
				bufferedGMLWriter.newLine();
				bufferedGMLWriter.write("\t\tid " + Math.abs(i));
				bufferedGMLWriter.newLine();
				bufferedGMLWriter.write("\t]");
				bufferedGMLWriter.newLine();
			}
			HashSet<Integer> parents = adjList.getParents(i);
			HashSet<Integer> children = adjList.getChildren(i);
			for (Integer parent : parents) {
				if (!written.contains(Math.abs(parent))) {
					written.add(Math.abs(parent));
					bufferedGMLWriter.write("\tnode [");
					bufferedGMLWriter.newLine();
					bufferedGMLWriter.write("\t\tid " + Math.abs(parent));
					bufferedGMLWriter.newLine();
					bufferedGMLWriter.write("\t]");
					bufferedGMLWriter.newLine();
				}
				sb.append("\tedge [");
				sb.append(System.getProperty("line.separator"));
				sb.append("\t\tsource ");
				sb.append(Math.abs(parent));
				sb.append(System.getProperty("line.separator"));
				sb.append("\t\ttarget ");
				sb.append(Math.abs(i));
				sb.append(System.getProperty("line.separator"));
				if (parent >= 0) {
					sb.append("\t\trepression ");
					sb.append(KeyParam.CYTOSCAPE_FALSE);
					sb.append(System.getProperty("line.separator"));
				} else {
					sb.append("\t\trepression ");
					sb.append(KeyParam.CYTOSCAPE_TRUE);
					sb.append(System.getProperty("line.separator"));
				}
				sb.append("\t]");
				sb.append(System.getProperty("line.separator"));
			}
			for (Integer child : children) {
				if (!written.contains(Math.abs(child))) {
					written.add(Math.abs(child));
					bufferedGMLWriter.write("\tnode [");
					bufferedGMLWriter.newLine();
					bufferedGMLWriter.write("\t\tid " + Math.abs(child));
					bufferedGMLWriter.newLine();
					bufferedGMLWriter.write("\t]");
					bufferedGMLWriter.newLine();
				}
				sb.append("\tedge [");
				sb.append(System.getProperty("line.separator"));
				sb.append("\t\tsource ");
				sb.append(Math.abs(i));
				sb.append(System.getProperty("line.separator"));
				sb.append("\t\ttarget ");
				sb.append(Math.abs(child));
				sb.append(System.getProperty("line.separator"));
				sb.append("\t\trepression ");
				sb.append(KeyParam.CYTOSCAPE_FALSE);
				sb.append(System.getProperty("line.separator"));
				sb.append("\t]");
				sb.append(System.getProperty("line.separator"));
			}
		}
		bufferedGMLWriter.write(sb.toString());
		bufferedGMLWriter.write("]");
		bufferedGMLWriter.flush();
		bufferedGMLWriter.close();
	}

	private static void writeGMLPreamble(BufferedWriter bufferedGMLWriter, ProblemInstance instance, String fileName,
										 String nameModifierString) throws IOException {
		bufferedGMLWriter.write("graph [");
		bufferedGMLWriter.newLine();
		bufferedGMLWriter.write("\tcomment \"" + "\"");
		bufferedGMLWriter.newLine();
		bufferedGMLWriter.write("\tdirected 1");
		bufferedGMLWriter.newLine();
	}

	public static void writeFullGraph(ArrayList<Future<double[]>> results, ProblemInstance instance) {
		String fileName = new SimpleDateFormat("yyyyMMddhhmm'.gml'").format(new Date());
		BufferedWriter bufferedGMLWriter = null;
		BufferedWriter bufferedNodeAttributeWriter = null;
		BufferedWriter bufferedPathwayWriter = null;
		BufferedWriter bufferedProteinWriter = null;
		boolean activesonly = true;
		try {
			bufferedNodeAttributeWriter = new BufferedWriter(new FileWriter(new File(instance.getOutputDirectory()
					.getPath() +
					System.getProperty("file.separator") + "NodeAttributes" + fileName)));
			HashSet<Integer> interactionSet = new HashSet<Integer>();
			HashMap<Integer, Object> integerObjectHashMap = instance.getIntegerObjectHashMap();
			for (Integer entityInteger : integerObjectHashMap.keySet()) {
				if (integerObjectHashMap.get(entityInteger) instanceof Conversion) {
					if (activesonly) {
						if (results.get(0).get()[entityInteger] != 0.0) {
							interactionSet.add(entityInteger);
						}
					} else {
						interactionSet.add(entityInteger);
					}
				}
			}
			bufferedGMLWriter = new BufferedWriter(new FileWriter(new File(instance.getOutputDirectory().getPath() +
					System.getProperty("file.separator") + "Network_" + "combined" + "_" + fileName)));
			writeGMLPreamble(bufferedGMLWriter, instance, fileName, "combined");
			writeGMLFile(bufferedGMLWriter, instance, interactionSet, fileName, "combined");
			writePathwayFile(bufferedPathwayWriter, instance, interactionSet, fileName, "combined", instance.getOutputDirectory());
			writeProteinList(bufferedProteinWriter, instance, interactionSet, fileName, "combined", instance.getOutputDirectory());
			writeFullNodeAttributeFile(bufferedNodeAttributeWriter, instance, interactionSet, results);
		} catch (Exception e) {
			log.error("ERROR : Differential writer error", e);
		} finally {
			IOUtils.closeQuietly(bufferedGMLWriter);
			IOUtils.closeQuietly(bufferedNodeAttributeWriter);
			IOUtils.closeQuietly(bufferedPathwayWriter);
			IOUtils.closeQuietly(bufferedProteinWriter);
		}
	}

	private static void writeFullNodeAttributeFile(BufferedWriter bufferedNodeAttributeWriter, ProblemInstance
			instance,
												   HashSet<Integer> interactionSet, ArrayList<Future<double[]>>
			results) throws Exception {
		HashSet<Integer> linesToWrite = new HashSet<Integer>();
		AdjList adjList = instance.getNetworkAdjList();
		for (Integer i : interactionSet) {
			linesToWrite.add(Math.abs(i));
			linesToWrite.addAll(adjList.getNeighbours(Math.abs(i)));
		}
		ArrayList<HashSet<Integer>> interactionLocationList = upstreamReactionTraverser(instance, results.get(0).get());
		HashMap<Integer, Object> integerObjectHashMap = instance.getIntegerObjectHashMap();
		bufferedNodeAttributeWriter.write("ID\tDISPLAY_NAME\tRDF_ID\tTYPE\tCELLULAR_LOCATION\tURL\tACTIVE" +
				"\tLOCATION" +
				"\tSETTER"+KeyParam.NEWLINE);
		for (Integer currentInteger : linesToWrite) {
			StringBuilder stringBuilder = new StringBuilder();
			Entity currentEntity = (Entity) integerObjectHashMap.get(Math.abs(currentInteger));
			stringBuilder.append(currentInteger);
			stringBuilder.append("\t");
			String name = getNodeName(currentEntity);
			stringBuilder.append(name);
			stringBuilder.append("\t");
			stringBuilder.append(currentEntity.getRDFId().replaceAll(instance.getModel().getXmlBase(), ""));
			stringBuilder.append("\t");
			if (currentEntity instanceof Interaction) {
				stringBuilder.append("INTERACTION\t");
			} else if (currentEntity instanceof Protein) {
				stringBuilder.append("PROTEIN\t");
			} else if (currentEntity instanceof Complex) {
				stringBuilder.append("COMPLEX\t");
			} else if (currentEntity instanceof Control) {
				stringBuilder.append("CONTROL\t");
			} else if (currentEntity instanceof SmallMolecule) {
				stringBuilder.append("SMALL_MOLECULE\t");
			} else if (currentEntity instanceof PhysicalEntity) {
				stringBuilder.append("PHYSICAL_ENTITY\t");
			} else {
				stringBuilder.append("OUTPUT_ERROR_UNCAUGHT_TYPE\t");
			}
			//Add cellular location if PhysicalEntity, otherwise skip (Interactions don't have locations):
			if (currentEntity instanceof Interaction) {
				stringBuilder.append("\t");
			} else {
				CellularLocationVocabulary cellularLocationVocabulary = ((PhysicalEntity) currentEntity)
						.getCellularLocation();
				if (cellularLocationVocabulary != null) {
					String locationString = cellularLocationVocabulary.getTerm().iterator().next();
					if (locationString == null) {
						stringBuilder.append(cellularLocationVocabulary.toString());
					} else {
						stringBuilder.append(locationString);
					}
				} else {
					//Logger warns in getNodeName method if there's no cellular location associated.
					stringBuilder.append("Not Recorded");
				}
				stringBuilder.append("\t");
			}
			//Add URL link to DB to node:
			if (currentEntity instanceof Interaction) {
				boolean appended = false;
				if (currentEntity.getXref() != null) {
					HashSet<Xref> xrefs = new HashSet<Xref>(currentEntity.getXref());
					for (Xref xr : xrefs) {
						if (xr != null && xr.getDb() != null && xr.getDb().equals("Reactome")) {
							stringBuilder.append(KeyParam.REACT_URL);
							stringBuilder.append(xr.getId());
							appended=true;
						}
					}
				}
				stringBuilder.append("\t");

			} else {
				stringBuilder.append("\t");
			}
			//Write the condition states for the node:
			for (int index = 0; index < results.size(); index++) {
				if (results.get(index).get()[Math.abs(currentInteger)] > KeyParam.SOLVER_INACTIVE) {
					stringBuilder.append(true);
				} else {
					stringBuilder.append(false);
				}
				stringBuilder.append("\t");
			}
			if (interactionLocationList.get(0).contains(currentInteger)) {
				stringBuilder.append("INTERNAL");
			} else if (interactionLocationList.get(1).contains(currentInteger)) {
				stringBuilder.append("UPSTREAM");
			} else if (interactionLocationList.get(2).contains(currentInteger)) {
				stringBuilder.append("DOWNSTREAM");
			}
			stringBuilder.append("\t");
			if (!(currentEntity instanceof Interaction) && instance.getSetterMap().containsKey(currentEntity)) {
				int settercode = instance.getSetterMap().get(currentEntity);
				if (settercode == KeyParam.SETTER_USER_SPEC) {
					stringBuilder.append("USER");
				} else if (settercode == KeyParam.SETTER_COMPLEX_FILL) {
					stringBuilder.append("COMPLEX FILL");
				} else if (settercode == KeyParam.SETTER_ESET_FILL) {
					stringBuilder.append("ESET FILL");
				}
			}
			stringBuilder.append(System.getProperty("line.separator"));
			bufferedNodeAttributeWriter.write(stringBuilder.toString());
		}
		bufferedNodeAttributeWriter.flush();
		bufferedNodeAttributeWriter.close();
	}

	private static void writeAllDownstream(BufferedWriter downstreamWriter, BufferedWriter downstreamProteinWriter,
										   ProblemInstance instance,
										   HashSet<Integer> interactionSet) throws IOException {
		LinkedBlockingQueue<Integer> processQueue = new LinkedBlockingQueue<Integer>();
		HashSet<PhysicalEntity> toWrite = new HashSet<PhysicalEntity>();
		HashMap<Integer, Object> integerObjectHashMap = instance.getIntegerObjectHashMap();
		HashMap<Object, Integer> objectIntegerHashMap = instance.getObjectIntegerHashMap();
		processQueue.addAll(interactionSet);
		HashSet<Integer> done = new HashSet<Integer>();
		while (!processQueue.isEmpty()) {
			Integer current = Math.abs(processQueue.poll());
			Object currentObj = integerObjectHashMap.get(current);
			if (currentObj instanceof PhysicalEntity) {
				PhysicalEntity pe = (PhysicalEntity) currentObj;
				for (Interaction i : pe.getParticipantOf()) {
					if (i instanceof Conversion && objectIntegerHashMap.containsKey(i) && !done.contains
							(objectIntegerHashMap.get(i))) {
						processQueue.add(objectIntegerHashMap.get(i));
						done.add(objectIntegerHashMap.get(i));
					}
				}
				for (Control c : pe.getControllerOf()) {
					for (Process p : c.getControlled()) {
						if (p instanceof Conversion && objectIntegerHashMap.containsKey(p) && !done.contains
								(objectIntegerHashMap.get(p))) {
							processQueue.add(objectIntegerHashMap.get(p));
							done.add(objectIntegerHashMap.get(p));
						}
					}
				}
			} else {
				Conversion conv = (Conversion) currentObj;
				for (PhysicalEntity rightPE : conv.getRight()) {
					if (!done.contains(objectIntegerHashMap.get(rightPE))) {
						processQueue.add(objectIntegerHashMap.get(rightPE));
						toWrite.add(rightPE);
					}
				}
			}
		}
		for (PhysicalEntity pe : toWrite) {
			downstreamWriter.write(pe.getDisplayName());
			downstreamWriter.write("\t");
			downstreamWriter.write(pe.getRDFId());
			downstreamWriter.newLine();
		}
		LinkedBlockingQueue<PhysicalEntity> entityLinkedBlockingQueue = new LinkedBlockingQueue<PhysicalEntity>();
		HashSet<Entity> doneProteinWriting = new HashSet<Entity>();
		entityLinkedBlockingQueue.addAll(toWrite);
		while (!entityLinkedBlockingQueue.isEmpty()) {
			PhysicalEntity entity = entityLinkedBlockingQueue.poll();
			if (!doneProteinWriting.contains(entity)) {
				if (entity instanceof Protein) {
					entityLinkedBlockingQueue.addAll(entity.getMemberPhysicalEntity());
					Protein protein = (Protein) entity;
					downstreamProteinWriter.write(getNodeName(protein));
					downstreamProteinWriter.write("\t");
					downstreamProteinWriter.write(Util.getUniProt(protein));
					downstreamProteinWriter.write("\t");
					downstreamProteinWriter.newLine();
				} else if (entity instanceof PhysicalEntity && !entity.getMemberPhysicalEntity().isEmpty()) {
					entityLinkedBlockingQueue.addAll(entity.getMemberPhysicalEntity());
				} else if (entity instanceof Complex) {
					entityLinkedBlockingQueue.addAll(((Complex) entity).getComponent());
				}
				doneProteinWriting.add(entity);
			}
		}
		downstreamProteinWriter.flush();
		downstreamWriter.flush();
	}

	public static HashSet<Integer> generatePredictionOutputSets(HashSet<Integer> highQualityPredictionSet,
													ProblemInstance instance, boolean limited) {
		HashSet<Integer> returnSet = new HashSet<Integer>();
		if (limited) {
			int currentdepth = KeyParam.TRAVERSAL_MED_CONF;
			for (Integer integer : highQualityPredictionSet) {
				returnSet.addAll(depthLimitedTraversal(integer, instance.getNetworkAdjList(), currentdepth, KeyParam
						.TRAVERSAL_MED_CONF));
			}
		} else {
			HashSet<Integer> fullTraverse = fullTraverse(highQualityPredictionSet, instance.getNetworkAdjList());
			HashMap<Integer, Object> ioMap = instance.getIntegerObjectHashMap();
			for (Integer idx : fullTraverse) {
				if (ioMap.get(idx) instanceof Interaction) returnSet.add(idx);
			}
		}
		return returnSet;
	}

	private static HashSet<Integer> depthLimitedTraversal(int i, AdjList adjList, int depth, int startdepth) {
		HashSet<Integer> newSet = new HashSet<Integer>();
		if ((startdepth-depth)%2==0) {
			newSet.add(i);
		}
		if (depth==0) {
			return newSet;
		} else if (depth > 0) {
			for (Integer childidx : adjList.getChildren(i)) {
				newSet.addAll(depthLimitedTraversal(Math.abs(childidx), adjList, depth - 1, startdepth));
			}
		}
		return newSet;
	}

	private static HashSet<Integer> fullTraverse(HashSet<Integer> integerSet, AdjList adjList) {
		try {
			LinkedBlockingQueue<Integer> processQueue = new LinkedBlockingQueue<Integer>();
			processQueue.addAll(integerSet);
			HashSet<Integer> done = new HashSet<Integer>();
			done.addAll(integerSet);
			while (!processQueue.isEmpty()) {
				Integer current = Math.abs(processQueue.poll());
				for (Integer childInteger : adjList.getChildren(current)) {
					if (!done.contains(childInteger)) {
						processQueue.add(childInteger);
						done.add(childInteger);
					}
				}
			}
			return done;
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
		//this shouldn't be reachable - compiler complains though.
		System.err.println("Unexpected termination in NetworkOutput.fullTraverse()");
		System.exit(-1);
		return new HashSet<Integer>();
	}

}
