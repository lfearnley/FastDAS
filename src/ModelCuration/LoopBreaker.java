package ModelCuration;

import FastDAS.KeyParam;
import org.apache.log4j.Logger;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.*;
import org.biopax.paxtools.model.level3.Process;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Stack;

/**
 * Created with IntelliJ IDEA.
 * User: uqlfearn
 * Date: 30/09/13
 * Time: 9:52 AM
 * To change this template use File | Settings | File Templates.
 */
public class LoopBreaker {

	private static org.apache.log4j.Logger log = Logger.getLogger(LoopBreaker.class);

	public static void breakLoops(Model model, File outputDirectory) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outputDirectory.getAbsolutePath() + System
					.getProperty("file.separator") + "loopbreakerLog.txt")));
			ArrayList<HashSet<BioPAXElement>> SCCs = tarjan(model);
			for (HashSet<BioPAXElement> SCC : SCCs) {
				processSCC(SCC, bw);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static void processSCC(HashSet<BioPAXElement> scc, BufferedWriter bufferedWriter) throws IOException {
		ArrayList<BioPAXElement> currentSCC = null;
		while (!(currentSCC = findLoopInSCC(scc)).isEmpty()) {
			bufferedWriter.newLine();
			bufferedWriter.write("------------------------------------------------------");
			bufferedWriter.newLine();
			//We're looking for the smallest element possible in the loop.
			//Most loops are caused by degradation or disassociation of complexes. As an example, the system:
			//A + B -> A:B
			//A:B + X -> A:B-p + Y
			//A:B-p -> A + B-p
			//The A:B complex can't form (and subsequently be modified) without the presence of A. In the A:B-p -> A +
			// B-p
			//reaction, the A:B-p -> A information is therefore redundant (it can't get more or less active in a
			// Boolean
			//formation, only be present or absent), and we want to remove it.
			//There's one of five cases:
			//1) X -> X - this generally is mislabelled/badly curated catalysis. Shouldn't happen, but we can test it
			// easily.
			//2) There's a control element in the SCC (ie feedback inhibition or upregulation) - we can't capture
			// inhibition loops
			//and we don't make things 'more on' in a Boolean system, so we break these.
			//3) {A:B:C:D} -> {A:B} - ie, the complex dissociates into other complexes (in which case this happens at
			// the smallest
			//complex.
			//4) {A:B:C:D} -> A - ie, the complex dissociates and spits a single protein out.
			//5) A-p -> A - ie, a modified protein is being stripped of its PTMs
			//Case 1: X->X (this is characterised by a SCC size of 2.
			if (currentSCC.size() == 2) {
				PhysicalEntity pe = currentSCC.get(0) instanceof PhysicalEntity ? (PhysicalEntity) currentSCC.get(0) :
						(PhysicalEntity) currentSCC.get(1);
				Conversion conv = currentSCC.get(0) instanceof Conversion ? (Conversion) currentSCC.get(0) :
						(Conversion) currentSCC.get(1);
				//Case 1A) A -> Rxn -> A : remove A from RHS of Rxn.
				if (conv.getLeft().contains(pe)) {
					conv.removeRight(pe);
					//Case 1B) A controls Rxn -> A : remove control of Rxn
				} else {
					HashSet<Control> ctrlset = new HashSet<Control>(pe.getControllerOf());
					for (Control c : ctrlset) {
						if (c.getControlled().contains(conv)) {
							c.removeControlled(conv);
						}
					}
				}
				continue;
			} else {
				//Initialise variables for iteration over the SCC.
				int startpoint = 0;
				int iteratorIndex = 0;
				BioPAXElement currentSmallest = null;
				BioPAXElement generator = null;
				BioPAXElement precedingElement = null;
				//Find the first PE and its generating
				if (currentSCC.get(0) instanceof Interaction) {
					startpoint = 1;
					generator = currentSCC.get(0);
					currentSmallest = currentSCC.get(1);
					precedingElement = currentSCC.get(currentSCC.size() - 1);
					iteratorIndex = 3;
				} else if (currentSCC.get(0) instanceof PhysicalEntity) {
					startpoint = 0;
					generator = currentSCC.get(currentSCC.size() - 1);
					currentSmallest = currentSCC.get(0);
					precedingElement = currentSCC.get(currentSCC.size() - 2);
					iteratorIndex = 2;
				}
				//Case 2: Control loop
				boolean iscontrolled = false;
				for (int i = startpoint; i < currentSCC.size(); i += 2) {
					PhysicalEntity pe = (PhysicalEntity) currentSCC.get(i);
					Conversion forward = i == (currentSCC.size() - 1) ? (Conversion) currentSCC.get(0) : (Conversion)
							currentSCC.get(i + 1);
					//If the left hand side of the next interaction doesn't contain the current element, it's a
					// controller.
					if (!forward.getLeft().contains(currentSCC.get(i))) {
						//Find the control element that's controlling this
						HashSet<Control> controlSet = new HashSet<Control>(pe.getControllerOf());
						for (Control ctrl : controlSet) {
							if (ctrl.getControlled().contains(forward)) {
								ctrl.removeControlled(forward);
								iscontrolled = true;
							}
						}
					}
				}
				if (iscontrolled) {
					continue;
				}
				//Case 3/4/5: {A:B:C:D} -> {A:B}||{A:B:C:D} -> A
				while (iteratorIndex < currentSCC.size()) {
					bufferedWriter.write("SCC IS : " + currentSCC.size());
					bufferedWriter.newLine();
					bufferedWriter.write("ITERATOR INDEX IS : " + iteratorIndex);
					bufferedWriter.newLine();
					PhysicalEntity currentElement = (PhysicalEntity) currentSCC.get(iteratorIndex);
					bufferedWriter.write("ELEMENT AT ITERATORINDEX IS : " + currentSCC.get(iteratorIndex).getRDFId());
					bufferedWriter.newLine();
					BioPAXElement bpe = iteratorIndex - 1 >= 0 ? currentSCC.get((iteratorIndex - 1)) : currentSCC.get(
							(currentSCC.size() + (iteratorIndex - 1)));
					bufferedWriter.write("GENERATOR IS : " + bpe.getRDFId());
					bufferedWriter.newLine();
					Conversion currentGenerator = iteratorIndex - 1 >= 0 ? (Conversion) currentSCC.get((iteratorIndex
							- 1)) : (Conversion) currentSCC.get((currentSCC.size() + (iteratorIndex - 1)));
					PhysicalEntity currentPrecedingElement = iteratorIndex - 2 >= 0 ? (PhysicalEntity) currentSCC.get(
							(iteratorIndex - 2)) : (PhysicalEntity) currentSCC.get((currentSCC.size() + (iteratorIndex
							- 2)));
					iteratorIndex += 2;
					if (currentElement instanceof Complex && currentPrecedingElement instanceof Complex && ((Complex)
							currentElement).getComponent().size() < ((Complex) currentPrecedingElement).getComponent()
							.size()) {
						//Case 3: Found a Complex -> smaller Complex dissociation event.
						if (currentSmallest instanceof Complex && ((Complex) currentElement).getComponent().size() < (
								(Complex) currentSmallest).getComponent().size()) {
							currentSmallest = currentElement;
							generator = currentGenerator;
							precedingElement = currentPrecedingElement;
						}
					} else if (currentElement instanceof Protein && currentPrecedingElement instanceof Complex) {
						//Case 4: Found a complex->protein dissociation event.
						currentSmallest = currentElement;
						generator = currentGenerator;
						precedingElement = currentPrecedingElement;
						break;
					} else if (currentElement instanceof Protein && currentPrecedingElement instanceof Protein &&
							currentElement.getFeature().size() < currentPrecedingElement.getFeature().size()) {
						//Case 5: Found a PTM-stripping event.
						if (currentSmallest instanceof Protein && currentElement.getFeature().size() < ((Protein)
								currentSmallest).getFeature().size()) {
							currentSmallest = currentElement;
							generator = currentGenerator;
							precedingElement = currentPrecedingElement;
						}
					}
				}
				//Break at generator->currentSmallest
				Conversion generatingConversion = (Conversion) generator;
				bufferedWriter.write(generatingConversion.getRDFId() + " removing " + currentSmallest.getRDFId() + " "
						+ generatingConversion.getRight().contains(currentSmallest));
				bufferedWriter.newLine();
				generatingConversion.removeRight((PhysicalEntity) currentSmallest);
				bufferedWriter.write("SCC WAS: ");
				bufferedWriter.newLine();
				for (BioPAXElement element : currentSCC) {
					if (element instanceof PhysicalEntity) {
						bufferedWriter.write(element.getRDFId() + "\t" + ((Entity) element).getDisplayName());
						bufferedWriter.newLine();
					} else if (element instanceof Conversion) {
						Conversion elConv = (Conversion) element;
						bufferedWriter.write(elConv.getRDFId() + "\t");
						for (PhysicalEntity pe : elConv.getLeft()) {
							bufferedWriter.write(pe.getDisplayName() + " ");
						}
						bufferedWriter.write(elConv.getDisplayName() + " ");
						for (PhysicalEntity pe : elConv.getRight()) {
							bufferedWriter.write(pe.getDisplayName() + " ");
						}
						bufferedWriter.newLine();
					}
				}
			}
		}
	}

	private static ArrayList<BioPAXElement> findLoopInSCC(HashSet<BioPAXElement> sccSet) {
		ArrayList<BioPAXElement> result = new ArrayList<BioPAXElement>();
		Stack<BioPAXElement> stack = new Stack<BioPAXElement>();
		stack.push(sccSet.iterator().next());
		HashMap<BioPAXElement, Integer> statusMap = new HashMap<BioPAXElement, Integer>();
		HashMap<BioPAXElement, BioPAXElement> parentMap = new HashMap<BioPAXElement, BioPAXElement>();
		while (!stack.isEmpty()) {
			BioPAXElement current = stack.pop();
			//If the node hasn't been visited:
			if (!statusMap.containsKey(current)) {
				//Set the node as visited, not closed
				statusMap.put(current, KeyParam.OPEN);
				stack.push(current);
				HashSet<BioPAXElement> successorElements = new HashSet<BioPAXElement>();
				if (current instanceof Conversion) {
					Conversion conv = (Conversion) current;
					for (PhysicalEntity pe : conv.getRight()) {
						if (sccSet.contains(pe)) {
							successorElements.add(pe);
						}
					}
				} else {
					PhysicalEntity pe = (PhysicalEntity) current;
					for (Interaction i : pe.getParticipantOf()) {
						if (i instanceof Conversion && !(i instanceof Degradation)) {
							Conversion conv = (Conversion) i;
							if (conv.getLeft().contains(pe) && sccSet.contains(conv)) {
								successorElements.add(conv);
							}
						}
					}
					for (Control ctrl : pe.getControllerOf()) {
						for (Process p : ctrl.getControlled()) {
							if (sccSet.contains(p)) {
								successorElements.add(p);
							}
						}
					}
				}
				for (BioPAXElement successor : successorElements) {
					if (!statusMap.containsKey(successor)) {
						parentMap.put(successor, current);
						stack.push(successor);
					} else if (statusMap.get(successor) == KeyParam.OPEN) {
						//CYCLE DETECTED:
						BioPAXElement pointer = current;
						while (!result.contains(successor)) {
							result.add(pointer);
							pointer = parentMap.get(pointer);
						}
						ArrayList<BioPAXElement> returnArrayList = new ArrayList<BioPAXElement>();
						for (int i = result.size() - 1; i >= 0; i--) {
							returnArrayList.add(result.get(i));
						}
						return returnArrayList;
					}
				}
			} else if (statusMap.get(current) == KeyParam.OPEN) {
				statusMap.put(current, KeyParam.CLOSED);
			}
		}
		return new ArrayList<BioPAXElement>();
	}


	/**
	 * Adapted from
	 * http://algowiki.net/wiki/index.php?title=Tarjan%27s_algorithm for use
	 * with integer adj matrix.
	 *
	 * @return
	 */
	private static ArrayList<HashSet<BioPAXElement>> tarjan(Model model) {
		int tarjIdx = 0;
		ArrayList<BioPAXElement> tarjStack = new ArrayList<BioPAXElement>();
		ArrayList<HashSet<BioPAXElement>> SCC = new ArrayList<HashSet<BioPAXElement>>();
		HashMap<BioPAXElement, Integer> tarjLowLinks = new HashMap<BioPAXElement, Integer>();
		HashMap<BioPAXElement, Integer> tarjIndexes = new HashMap<BioPAXElement, Integer>();
		//Initialise variables in HashMaps:
		for (PhysicalEntity pe : model.getObjects(PhysicalEntity.class)) {
			tarjLowLinks.put(pe, -1);
			tarjIndexes.put(pe, -1);
		}
		for (Conversion conv : model.getObjects(Conversion.class)) {
			if (!(conv instanceof Degradation)) {
				tarjLowLinks.put(conv, -1);
				tarjIndexes.put(conv, -1);
			}
		}
		tarjIdx = 0;
		for (PhysicalEntity pe : model.getObjects(PhysicalEntity.class)) {
			if (tarjIndexes.get(pe) == -1) {
				tarjanSub(pe, model, tarjIndexes, tarjLowLinks, tarjStack, SCC, tarjIdx);
			}
		}
		return SCC;
	}

	private static void tarjanSub(BioPAXElement element, Model model, HashMap<BioPAXElement, Integer> tarjIndexes,
								  HashMap<BioPAXElement, Integer> tarjLowLinks, ArrayList<BioPAXElement> tarjStack,
								  ArrayList<HashSet<BioPAXElement>> SCC, int tarjIdx) {
		tarjIndexes.put(element, tarjIdx);
		tarjLowLinks.put(element, tarjIdx);
		tarjIdx++;
		tarjStack.add(0, element);
		HashSet<BioPAXElement> rightNeighbours = new HashSet<BioPAXElement>();
		if (element instanceof PhysicalEntity) {
			PhysicalEntity pe = (PhysicalEntity) element;
			for (Interaction i : pe.getParticipantOf()) {
				if (i instanceof Conversion && !(i instanceof Degradation)) {
					Conversion conv = (Conversion) i;
					if (conv == null) {
						continue;
					}
					if (conv.getLeft().contains(pe) && tarjIndexes.containsKey(conv)) {
						rightNeighbours.add(conv);
					}
				}
			}
			for (Control ctrl : pe.getControllerOf()) {
				if (tarjIndexes.containsKey(ctrl)) {
					if (ctrl == null) {
						continue;
					}
					for (Process p : ctrl.getControlled()) {
						if (tarjIndexes.containsKey(p)) {
							rightNeighbours.add(p);
						}
					}
				}
			}
		} else if (element instanceof Conversion) {
			Conversion conv = (Conversion) element;
			for (PhysicalEntity pe : conv.getRight()) {
				if (tarjIndexes.containsKey(pe)) {
					if (pe == null) {
						continue;
					}
					rightNeighbours.add(pe);
				}
			}
		}
		for (BioPAXElement nextElement : rightNeighbours) {
			if (tarjIndexes.get(nextElement) == -1) {
				tarjanSub(nextElement, model, tarjIndexes, tarjLowLinks, tarjStack, SCC, tarjIdx);
				tarjLowLinks.put(element, Math.min(tarjLowLinks.get(element), tarjLowLinks.get(nextElement)));
			} else if (tarjStack.contains(nextElement)) {
				tarjLowLinks.put(element, Math.min(tarjLowLinks.get(element), tarjIndexes.get(nextElement)));
			}
		}
		if (tarjLowLinks.get(element) == tarjIndexes.get(element)) {
			HashSet<BioPAXElement> singleSCC = new HashSet<BioPAXElement>();
			BioPAXElement node;
			do {
				node = tarjStack.remove(0);
				singleSCC.add(node);
			} while (!node.equals(element));
			if (singleSCC.size() > 1) {
				SCC.add(singleSCC);
			}
		}
	}

}
