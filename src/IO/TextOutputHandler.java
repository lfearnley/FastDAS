package IO;

import FastDAS.KeyParam;
import LogicalSystem.ProblemInstance;
import LogicalSystem.Tuple2;
import LogicalSystem.Util;
import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.*;
import org.biopax.paxtools.model.level3.Process;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.security.Key;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.regex.Matcher;

/**
 * Created with IntelliJ IDEA.
 * User: uqlfearn
 * Date: 16/08/13
 * Time: 1:04 PM
 * To change this template use File | Settings | File Templates.
 */
public class TextOutputHandler {

	private static org.apache.log4j.Logger log = Logger.getLogger(TextOutputHandler.class);

	public static void writePathwayHierarchy(Model model, File outputDirectory) {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(outputDirectory.getPath() + System.getProperty
					("file.separator") + "pathway.txt"));
			LinkedBlockingQueue<Pathway> pathwayQueue = new LinkedBlockingQueue<Pathway>();
			HashMap<Pathway, Boolean> pathwayWriteFlags = new HashMap<Pathway, Boolean>();
			for (Pathway path : model.getObjects(Pathway.class)) {
				if (path.getPathwayComponentOf().isEmpty()) {
					pathwayQueue.add(path);
				}
				pathwayWriteFlags.put(path, false);
			}
			while (!pathwayQueue.isEmpty()) {
				String tabString = "";
				Pathway current = pathwayQueue.poll();
				HashMap<Pathway, Boolean> localPathwayWriteFlags = new HashMap<Pathway, Boolean>();
				writePathwayWithIndent(current, current, pathwayWriteFlags, localPathwayWriteFlags, writer, tabString);
			}
			for (Pathway p : pathwayWriteFlags.keySet()) {
				if (!pathwayWriteFlags.get(p)) {
					System.err.println("Did not write : " + p.getDisplayName());
				}
			}
			writer.flush();
			writer.close();
		} catch (Exception e) {
			log.error("Error in writing pathway hierarchy file : ", e);
			log.error("Shutting down.");
			System.exit(-1);
		}
	}

	private static void writePathwayWithIndent(Pathway parent, Pathway poll,
											   HashMap<Pathway, Boolean> globalPathwayWriteFlags,
											   HashMap<Pathway, Boolean> localPathwayWriteFlags,
											   BufferedWriter writer, String indent) throws IOException,
			InterruptedException {
		String stableID = "unknown";
		if (poll.getXref() != null) {
			HashSet xrefs = new HashSet<Xref>(poll.getXref());
			if (!xrefs.isEmpty()) {
				for (Xref xref : poll.getXref()) {
					if (xref.getDb() == null) {
						continue;
					} else if (xref.getDb().equals("Reactome")) {
						stableID = xref.getId();
					}
				}
			}
		} else {
			//System.err.println(poll.getRDFId() + " has no Xrefs");
		}
		//writer.write(indent + stableID + "," + poll.getDisplayName() + "," + getPathwayComponentUniprots(poll));
		writer.write(indent + poll.getDisplayName());
		writer.newLine();
		globalPathwayWriteFlags.put(poll, true);
		localPathwayWriteFlags.put(poll, true);
		String newIndent = indent + "\t";
		for (Process p : poll.getPathwayComponent()) {
			if (p instanceof Pathway) {
				if (localPathwayWriteFlags.containsKey(p) && localPathwayWriteFlags.get(p)) {
					continue;
				} else {
					writePathwayWithIndent(parent, (Pathway) p, globalPathwayWriteFlags, localPathwayWriteFlags,
							writer, newIndent);
				}
			}
		}
	}

	private static String getPathwayComponentUniprots(Pathway pathway) throws InterruptedException {
		StringBuilder sb = new StringBuilder();
		HashSet<String> stringHashSet = new HashSet<String>();
		LinkedBlockingQueue<Entity> eQueue = new LinkedBlockingQueue<Entity>();
		eQueue.put(pathway);
		while (!eQueue.isEmpty()) {
			Entity polledEntity = eQueue.poll();
			if (polledEntity instanceof Pathway) {
				eQueue.addAll(((Pathway) polledEntity).getPathwayComponent());
			} else if (polledEntity instanceof Interaction) {
				eQueue.addAll(((Interaction) polledEntity).getParticipant());
			} else if (polledEntity instanceof Complex) {
				eQueue.addAll(((Complex) polledEntity).getMemberPhysicalEntity());
				eQueue.addAll(((Complex) polledEntity).getComponent());
			} else if (polledEntity instanceof Protein) {
				Protein current = (Protein) polledEntity;
				eQueue.addAll(current.getMemberPhysicalEntity());
				String uniProt = "";
				EntityReference entityReference = current.getEntityReference();
				if (entityReference != null) {
					//Check for uniprot in EntityReference:
					for (String s : entityReference.getName()) {
						Matcher matcher = KeyParam.UNIPROT_PATTERN_ONE.matcher(s);
						if (matcher.find()) {
							uniProt = matcher.group(1);
							stringHashSet.add(matcher.group(1));
							break;
						} else {
							matcher = KeyParam.UNIPROT_PATTERN_TWO.matcher(s);
							if (matcher.find()) {
								uniProt = matcher.group(1);
								stringHashSet.add(matcher.group(1));
								break;
							}
						}
					}
					if (uniProt.equals("")) {
						if (entityReference.getXref() != null) {
							HashSet<Xref> xrefs = new HashSet<Xref>(entityReference.getXref());
							for (Xref xref : xrefs) {
								if (xref.getDb().equals("UniProt")) {
									stringHashSet.add(xref.getId());
									break;
								}
							}
						} else {
							log.warn(current.getRDFId() + " has no Xrefs");
						}
					}
				} else {
					log.warn(current.getRDFId() + " has no EntityReference");
				}
			}
		}
		StringBuilder stringBuilder = new StringBuilder();
		for (String s : stringHashSet) {
			stringBuilder.append(s);
			stringBuilder.append(",");
		}
		return stringBuilder.toString();
	}

	public static void writeSaveFileForEditing(Model model, String outputString) {
		try {
			BufferedWriter speciesWriter = new BufferedWriter(new FileWriter(new File(outputString)));
			for (PhysicalEntity pe : model.getObjects(PhysicalEntity.class)) {
				//Skip species that don't participate in reactions - these will be removed.
				if (pe.getParticipantOf().isEmpty()) {
					continue;
				}
				boolean input = true;
				boolean output = true;
				for (Interaction i : pe.getParticipantOf()) {
					if (i instanceof Conversion) {
						Conversion conv = (Conversion) i;
						if (conv.getLeft().contains(pe)) {
							output = false;
						} else if (conv.getRight().contains(pe)) {
							input = false;
						}
					} else if (i instanceof Control) {
						output = false;
					}
				}
				StringBuilder sb = new StringBuilder();
				sb.append(pe.getClass());
				sb.append("\t");
				sb.append(pe.getRDFId().replaceAll(model.getXmlBase(), ""));
				sb.append("\t");
				sb.append(pe.getDisplayName());
				sb.append("\t");
				if (pe.getCellularLocation() != null) {
					sb.append(pe.getCellularLocation().toString());
				}
				sb.append("\t");
				sb.append(pe.getParticipantOf().size());
				sb.append("\t");
				if (input && output) {
					sb.append("SINGLETON");
				} else if (input) {
					sb.append("INPUT");
				} else if (output) {
					sb.append("OUTPUT");
				} else {
					sb.append("INTERNAL");
				}
				sb.append("\t");
				if (pe instanceof SmallMolecule) {
					sb.append("REMOVE");
				} else {
					sb.append("RETAIN");
				}
				sb.append(System.lineSeparator());
				speciesWriter.write(sb.toString());
			}
			speciesWriter.flush();
			speciesWriter.close();
			writeMappingFile(model);
		} catch (Exception e) {
			log.error("Error in writing species deletion file : ", e);
			log.error("Shutting down.");
			System.exit(-1);
		}

	}

	public static void writeMappingFile(Model model) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File("mapping.txt")));
			bw.write("RDF\tUniProt\tREACT_ID");
			bw.write(KeyParam.NEWLINE);
			for (PhysicalEntity pe : model.getObjects(PhysicalEntity.class)) {
				bw.write(pe.getRDFId().replaceAll(model.getXmlBase(),""));
				bw.write("\t");
				if (pe instanceof Protein) {
					bw.write(Util.getUniProt((Protein) pe));
				}
				bw.write("\t");
				if (pe.getXref() != null) {
					HashSet<Xref> xrefs = new HashSet<Xref>(pe.getXref());
					for (Xref xr : xrefs) {
						if (xr != null && xr.getDb() != null && xr.getDb().equals("Reactome")) {
							bw.write(xr.getId());
						}
					}
				}
				bw.write("\t");
				if (pe == null || pe.getDisplayName() == null) {
					bw.write("\t");
				} else {
					bw.write(pe.getDisplayName());
				}
				bw.write(KeyParam.NEWLINE);
			}
			bw.flush();
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void writeSIFFromNetwork(Model model, String filename) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(filename)));
			for (Conversion c : model.getObjects(Conversion.class)) {
				for (PhysicalEntity pe : c.getLeft()) {
					bw.write(pe.getRDFId() + "\tpp\t" + c.getRDFId());
					bw.newLine();
				}
				for (PhysicalEntity pe : c.getRight()) {
					bw.write(c.getRDFId() + "\tpp\t" + pe.getRDFId());
					bw.newLine();
				}
			}
			for (Control ctrl : model.getObjects(Control.class)) {
				if (ctrl.getControlType().equals(ControlType.INHIBITION) || ctrl.getControlType().equals(ControlType
						.INHIBITION_OTHER) ||
						ctrl.getControlType().equals(ControlType.INHIBITION_ALLOSTERIC) || ctrl.getControlType()
						.equals(ControlType.INHIBITION_COMPETITIVE)
						|| ctrl.getControlType().equals(ControlType.INHIBITION_IRREVERSIBLE) || ctrl.getControlType()
						.equals(ControlType.INHIBITION_NONCOMPETITIVE)
						|| ctrl.getControlType().equals(ControlType.INHIBITION_UNCOMPETITIVE)) {
					if (ctrl.getControlled().iterator().hasNext()) {
						bw.write(ctrl.getRDFId() + "\tinh\t" + ctrl.getControlled().iterator().next().getRDFId());
						bw.newLine();
					} else {
						//Control of removed element (will be stripped out).
						continue;
					}
				} else {
					if (ctrl.getControlled().iterator().hasNext()) {
						bw.write(ctrl.getRDFId() + "\tctrl\t" + ctrl.getControlled().iterator().next().getRDFId());
						bw.newLine();
					} else {
						//Control of removed element (will be stripped out).
						continue;
					}
				}
				for (Controller ctrlr : ctrl.getController()) {
					bw.write(ctrlr.getRDFId() + "\tctrlr\t" + ctrl.getRDFId());
					bw.newLine();
				}
			}
			bw.flush();
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void writeResultConditionComparison(ProblemInstance instance,
													  HashMap<Integer, ArrayList<Tuple2<HashSet<Integer>,
															  HashSet<Integer>>>> resultMap,
													  int comparisonIndex) {
		BufferedWriter bufferedResultWriter = null;
		try {
			bufferedResultWriter = new BufferedWriter(new FileWriter(instance.getOutputDirectory().getPath()
					+ System.getProperty("file.separator") + "comparisonOutput_Control-" + (comparisonIndex + 1) + "" +
					".txt"));
			for (Integer resultCount : resultMap.keySet()) {
				bufferedResultWriter.write(resultCount);
				bufferedResultWriter.write("\t");
				bufferedResultWriter.write(resultMap.get(resultCount).size());
				bufferedResultWriter.newLine();
			}
			bufferedResultWriter.flush();
		} catch (Exception e) {
			log.error("ERROR : Error writing comparison between Control and instance " + comparisonIndex, e);
		} finally {
			IOUtils.closeQuietly(bufferedResultWriter);
		}
		//To change body of created methods use File | Settings | File Templates.
	}

	public static void writeEffectCounts(LinkedBlockingQueue<ArrayList<Integer>> solsizes, File outputDirectory) {
		BufferedWriter outputWriter = null;
		try {
			outputWriter = new BufferedWriter(new FileWriter(new File(outputDirectory.getPath() + KeyParam
					.PATH_SEPARATOR + "productcounts.txt")));
			outputWriter.write("HIGH\tMED\tLOW");
			outputWriter.write(KeyParam.NEWLINE);
			while (!solsizes.isEmpty()) {
				ArrayList<Integer> iList = solsizes.poll();
				outputWriter.write(iList.get(0));
				outputWriter.write("\t");
				outputWriter.write(iList.get(1));
				outputWriter.write("\t");
				outputWriter.write(iList.get(2));
				outputWriter.write(KeyParam.NEWLINE);
			}
			outputWriter.flush();
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			IOUtils.closeQuietly(outputWriter);
		}
	}
}
