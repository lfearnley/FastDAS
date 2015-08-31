package LogicalSystem;

import FastDAS.KeyParam;
import IO.NetworkOutputHandler;
import org.apache.log4j.Logger;
import org.biopax.paxtools.impl.level3.PhysicalEntityImpl;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.*;
import org.biopax.paxtools.model.level3.Process;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.Future;
import java.util.regex.Matcher;

/**
 * Created with IntelliJ IDEA.
 * User: uqlfearn
 * Date: 23/08/13
 * Time: 2:26 PM
 * To change this template use File | Settings | File Templates.
 */
public class Util {

	private static org.apache.log4j.Logger log = Logger.getLogger(Util.class);


	public static void removeInteraction(Interaction interactionToRemove, Model m) {
		HashSet<Pathway> pathwaySet = new HashSet<Pathway>(interactionToRemove.getPathwayComponentOf());
		for (Pathway path : pathwaySet) {
			path.removePathwayComponent(interactionToRemove);
		}
		HashSet<Xref> xrefs = new HashSet<Xref>(interactionToRemove.getXref());
		for (Xref xref : xrefs) {
			m.remove(xref);
		}
		HashSet<Entity> entitySet = new HashSet<Entity>(interactionToRemove.getParticipant());
		for (Entity participant : entitySet) {
			interactionToRemove.removeParticipant(participant);
			if (participant instanceof PhysicalEntity) {
				PhysicalEntity physicalEntity = (PhysicalEntity) participant;
				if (physicalEntity.getComponentOf().isEmpty() && physicalEntity.getParticipantOf().isEmpty() &&
						physicalEntity.getMemberPhysicalEntityOf().isEmpty()) {
					m.remove(physicalEntity);
				}
			}
		}
		m.remove(interactionToRemove);
	}

	public static Control duplicateControl(Control control, Model m, long expansionCount, String xmlBase) {
		Control newControl = null;
		if (control instanceof Catalysis) {
			newControl = m.addNew(Catalysis.class, control.getRDFId().replaceAll(xmlBase, "") + "EXPAND_" +
					expansionCount);
			expansionCount++;
			Catalysis oldCatalysis = (Catalysis) control;
			Catalysis newCatalysis = (Catalysis) newControl;
			for (PhysicalEntity pe : oldCatalysis.getCofactor()) {
				newCatalysis.addCofactor(pe);
			}
			//Check to see if there are any non-co-factor controllers
			for (Controller pe : oldCatalysis.getController()) {
				newCatalysis.addController(pe);
			}
		} else {
			newControl = m.addNew(Control.class, control.getRDFId().replaceAll(xmlBase, "") + "EXPAND_" +
					expansionCount);
			for (Controller pe : control.getController()) {
				newControl.addController(pe);
			}
		}
		for (Process controlled : control.getControlled()) {
			newControl.addControlled(controlled);
		}
		for (Pathway pathway : control.getPathwayComponentOf()) {
			pathway.addPathwayComponent(newControl);
		}
		if (control.getControlType() != null) {
			newControl.setControlType(control.getControlType());
		}
		copyMetadata(control, newControl);
		return newControl;
	}

	private static void copyMetadata(Interaction original, Interaction copy) {
		if (original.getDisplayName() != null) {
			copy.setDisplayName(original.getDisplayName());
		}
		if (original.getStandardName() != null) {
			copy.setStandardName(original.getStandardName());
		}
		if (original.getXref() != null && !original.getXref().isEmpty()) {
			for (Xref xr : original.getXref()) {
				copy.addXref(xr);
			}
		}
		if (original.getName() != null && !original.getName().isEmpty()) {
			for (String s : original.getName()) {
				copy.addName(s);
			}
		}
	}

	private static void copyMetadata(PhysicalEntity original, PhysicalEntity copy) {
		if (original.getDisplayName() != null) {
			copy.setDisplayName(original.getDisplayName());
		}
		if (original.getStandardName() != null) {
			copy.setStandardName(original.getStandardName());
		}
		if (original.getXref() != null && !original.getXref().isEmpty()) {
			for (Xref xr : original.getXref()) {
				copy.addXref(xr);
			}
		}
		if (original.getName() != null && !original.getName().isEmpty()) {
			for (String s : original.getName()) {
				copy.addName(s);
			}
		}
	}

	public static void replaceInConversion(Conversion conversion, PhysicalEntity toreplace, PhysicalEntity
			replacement) {
		HashSet<PhysicalEntity> lhsPE = new HashSet<PhysicalEntity>(conversion.getLeft());
		for (PhysicalEntity pe : lhsPE) {
			if (pe.equals(toreplace)) {
				Stoichiometry participantStoichiometry = Util.getParticipantStoichiometry(conversion, toreplace);
				if (participantStoichiometry != null) {
					participantStoichiometry.setPhysicalEntity(replacement);
				}
				conversion.addLeft(replacement);
				conversion.removeLeft(toreplace);
			}
		}
		HashSet<PhysicalEntity> rhsPE = new HashSet<PhysicalEntity>(conversion.getRight());
		for (PhysicalEntity pe : rhsPE) {
			if (pe.equals(toreplace)) {
				Stoichiometry participantStoichiometry = Util.getParticipantStoichiometry(conversion, toreplace);
				if (participantStoichiometry != null) {
					participantStoichiometry.setPhysicalEntity(replacement);
				}
				conversion.addRight(replacement);
				conversion.removeRight(toreplace);
			}
		}
	}

	public static boolean hasComplexOverlap(Complex container, Complex contained) {
		if (contained.getComponent().isEmpty()) {
			return container.getComponent().contains(contained);
		}
		for (PhysicalEntity pe : contained.getComponent()) {
			if (!container.getComponent().contains(pe)) {
				for (PhysicalEntity containerPE : container.getComponent()) {
					boolean match = false;
					if (containerPE instanceof Protein && pe instanceof Protein) {
						Protein containerProt = (Protein) containerPE;
						Protein containedProt = (Protein) pe;
						if (containerProt.getEntityReference() != null && containedProt.getEntityReference() != null &&
								containerProt.getEntityReference().equals(containedProt.getEntityReference())) {
							match = true;
						}
					}
					if (!match) {
						return false;
					}
				}
			} else {
				Stoichiometry containerStoich = Util.getComponentStoichiometry(container, pe);
				Stoichiometry containedStoich = Util.getComponentStoichiometry(contained, pe);
				if (containerStoich != null && containedStoich != null) {
					if (containerStoich.getStoichiometricCoefficient() != containedStoich.getStoichiometricCoefficient
							()) {
						return false;
					}
				} else if ((containerStoich == null && containedStoich != null) || (containerStoich != null &&
						containedStoich == null)) {
					return false;
				} else {
					continue;
				}
			}
		}
		return true;
	}

	public static Conversion duplicateConversion(Conversion conversion, Model m, long expansionCount, String xmlBase) {
		Conversion newConversion = m.addNew(Conversion.class, conversion.getRDFId().replaceAll(xmlBase, "") +
				"EXPAND_" + expansionCount);
		expansionCount++;
		for (PhysicalEntity pe : conversion.getLeft()) {
			//Get corresponding stoichiometry:
			Stoichiometry participantStoichiometry = Util.getParticipantStoichiometry(conversion, pe);
			newConversion.addLeft(pe);
			if (participantStoichiometry != null) {
				newConversion.addParticipantStoichiometry(participantStoichiometry); //TODO check whether need to
				// duplicate stoich
			}

		}
		for (PhysicalEntity pe : conversion.getRight()) {
			//Get corresponding stoichiometry:
			Stoichiometry participantStoichiometry = Util.getParticipantStoichiometry(conversion, pe);
			newConversion.addRight(pe);
			if (participantStoichiometry != null) {
				newConversion.addParticipantStoichiometry(participantStoichiometry); //TODO check whether need to
				// duplicate stoich
			}
		}
		//Update pathways to contain copy.
		for (Pathway pathway : conversion.getPathwayComponentOf()) {
			pathway.addPathwayComponent(newConversion);
		}
		copyMetadata(conversion, newConversion);
		return newConversion;
	}

	public static String getUniProt(Protein protein) {
		String uniProt = "";
		EntityReference entityReference = protein.getEntityReference();
		if (entityReference != null) {
			//Check for uniprot in EntityReference:
			for (String s : entityReference.getName()) {
				Matcher matcher = KeyParam.UNIPROT_PATTERN_ONE.matcher(s);
				if (matcher.find()) {
					uniProt = matcher.group(1);
					break;
				}
			}
			if (uniProt.equals("")) {
				HashSet<Xref> xrefs = new HashSet<Xref>(entityReference.getXref());
				if (xrefs != null) {
					for (Xref xref : xrefs) {
						if (xref.getDb().equals("UniProt")) {
							uniProt = xref.getId();
							break;
						}
					}
					log.warn(protein.getRDFId() + " has no UniProts");
				} else {
					log.warn(protein.getRDFId() + " has no Xrefs");
				}
			}
		} else {
			log.warn(protein.getRDFId() + " has no EntityReference");
		}
		return uniProt;
	}

	/**
	 * Retrieves the Stoichiometry object associated with a specific component of a Complex.
	 *
	 * @param complex   A BioPAX Level 3 Complex object.
	 * @param component A PhysicalEntity contained within the component set of complex.
	 * @return The Stoichiometry associated with the component.
	 */
	public static Stoichiometry getComponentStoichiometry(Complex complex, PhysicalEntity component) {
		Stoichiometry stoic = null;
		try {
			Set<Stoichiometry> stoichSet = complex.getComponentStoichiometry();
			//Make sure there -is- a stoichiometry associated with the complex - assume 1 if not.
			if (stoichSet != null) {
				for (Stoichiometry componentStoichiometry : stoichSet) {
					if (componentStoichiometry.getPhysicalEntity().equals(component)) {
						stoic = componentStoichiometry;
						break;
					}
				}
			}
		} catch (Exception e) {
			stoic = null;
		} finally {
			return stoic;
		}
	}

	public static Stoichiometry getEquivalentStoichiometryFromSet(Stoichiometry stoichiometry, HashSet<Stoichiometry>
			stoichSet) {
		for (Stoichiometry candidateStoichiometry : stoichSet) {
			if (candidateStoichiometry.getStoichiometricCoefficient() == stoichiometry.getStoichiometricCoefficient
					() &&
					candidateStoichiometry.getPhysicalEntity().equals(stoichiometry.getPhysicalEntity())) {
				return candidateStoichiometry;
			}
		}
		return null;
	}

	/**
	 * Method to replace an entity with another in the model, deleting the entity once replacement takes place.
	 *
	 * @param entity      PhysicalEntity to be replaced.
	 * @param replacement Replacement PhysicalEntity
	 * @param model       Model containing the entity and its replacement
	 */
	public static void replaceEntity(PhysicalEntity entity, PhysicalEntity replacement, Model model) {
		HashSet<Interaction> interactionSet = new HashSet<Interaction>(entity.getParticipantOf());
		for (Interaction i : interactionSet) {
			if (i instanceof Control) {
				//Handle Catalysis differently due to cofactors (entity.getParticipantOf will send back cofactors as
				// well):
				if (i instanceof Catalysis) {
					Catalysis catalysis = (Catalysis) i;
					if (catalysis.getCofactor().contains(entity)) {
						catalysis.removeCofactor(entity);
						catalysis.addCofactor(replacement);
					}
					if (catalysis.getController().contains(entity)) {
						catalysis.removeController(entity);
						catalysis.addController(replacement);
					}
				} else {
					Control control = (Control) i;
					control.removeController(entity);
					control.addController(replacement);
				}
			} else if (i instanceof Conversion) {
				Conversion conversion = (Conversion) i;
				//If there's a reaction stoichiometry here then update that to the replacement entity
				Stoichiometry stoichiometry = Util.getParticipantStoichiometry(conversion, entity);
				if (stoichiometry != null) {
					stoichiometry.setPhysicalEntity(replacement);
				}
				//Check left hand side:
				if (conversion.getLeft().contains(entity)) {
					conversion.removeLeft(entity);
					conversion.addLeft(replacement);
				}
				//Check right hand side (note not 'else if' - in case of misrepresented catalysis putting entity on
				//both LHS and RHS, this preserves stoichiometric balance):
				if (conversion.getRight().contains(entity)) {
					conversion.removeRight(entity);
					conversion.addRight(replacement);
				}
			} else if (i instanceof MolecularInteraction) {
				log.warn("WARN : MolecularInteraction present " + i.getRDFId());
				System.exit(-1);
				//omitted TemplateReaction/GeneticInteraction because there should be no complexes involved in
				//these processes, only proteins/genes (respectively).
			}
		}
		//Replace the entity in any complexes it participates in.
		HashSet<Complex> complexes = new HashSet<Complex>(entity.getComponentOf());
		for (Complex c : complexes) {
			Stoichiometry componentStoich = Util.getComponentStoichiometry(c, entity);
			c.addComponent(replacement);
			componentStoich.setPhysicalEntity(replacement);
			c.removeComponent(entity);
		}
		//Replace the entity in any entitySets it participates in.
		HashSet<PhysicalEntity> entitySets = new HashSet<PhysicalEntity>(entity.getMemberPhysicalEntityOf());
		for (PhysicalEntity entitySet : entitySets) {
			entitySet.removeMemberPhysicalEntity(entity);
			entitySet.addMemberPhysicalEntity(replacement);
		}
		//All replacements are done, remove the entity from the model.
		model.remove(entity);
	}

	public static Stoichiometry getParticipantStoichiometry(Conversion conversion, PhysicalEntity entity) {
		Stoichiometry stoichiometry = null;
		try {
			for (Stoichiometry s : conversion.getParticipantStoichiometry()) {
				if (s.getPhysicalEntity().equals(entity)) {
					stoichiometry = s;
					break;
				}
			}
		} catch (Exception e) {
			stoichiometry = null;
		} finally {
			return stoichiometry;
		}
	}

	private static void replaceInInteraction(Interaction interaction, PhysicalEntity toReplace, PhysicalEntity
			replacement) {
		if (interaction instanceof Control) {
			Control control = (Control) interaction;
			if (control.getController().contains(toReplace)) {
				control.addController(replacement);
				control.removeController(toReplace);
			}
		} else if (interaction instanceof Conversion) {
			Conversion conversion = (Conversion) interaction;
			Stoichiometry stoich = null;
			for (Stoichiometry stoichiometry : conversion.getParticipantStoichiometry()) {
				if (stoichiometry.getPhysicalEntity().equals(toReplace)) {
					stoich = stoichiometry;
				}
			}
			if (conversion.getLeft().contains(interaction)) {
				conversion.addLeft(replacement);
				stoich.setPhysicalEntity(replacement);
				conversion.removeLeft(toReplace);
			}
			if (conversion.getRight().contains(interaction)) {
				conversion.addRight(replacement);
				stoich.setPhysicalEntity(replacement);
				conversion.removeRight(toReplace);
			}
		} else {
			log.error("Interaction of unknown type at replacement of an element : " + interaction.getRDFId());

		}
	}

	private PhysicalEntity getEquivalentFromSet(PhysicalEntity set, PhysicalEntity element) {
		//Case 1: set is single-level, contains exact element:
		if (set.getMemberPhysicalEntity().contains(element)) {
			return element;
			//Case 2: set contains element, but element has been modified:
		} else {
			//Iterate through each element in the set:
			for (PhysicalEntity pe : set.getMemberPhysicalEntity()) {
				if (pe instanceof Protein && isEquivalent(pe, element)) {
					return pe;
				}
			}
		}
		//
		return null;
	}

	private boolean isEquivalent(PhysicalEntity pe, PhysicalEntity element) {
		if (pe instanceof Protein && element instanceof Protein) {
			Protein candidateProtein = (Protein) pe;
			Protein elementProtein = (Protein) element;
			if (candidateProtein.getEntityReference() != null && elementProtein.getEntityReference() != null) {
				if (candidateProtein.getEntityReference().equals(elementProtein.getEntityReference())) {
					return true;
				}
			}
		}
		return false;
	}

	public static Complex getEquivalentComplex(Model model, Complex originalComplex, PhysicalEntity
			componentToReplace, PhysicalEntity replacement, long splitcount) {
		HashSet<BioPAXElement> newElements = new HashSet<BioPAXElement>();
		Complex newComplex = model.addNew(Complex.class, originalComplex.getRDFId().replaceAll(model.getXmlBase(), "")
				+ "_" + splitcount);
		newElements.add(newComplex);
		splitcount++;
		for (PhysicalEntity component : originalComplex.getComponent()) {
			Stoichiometry componentStoichiometry = Util.getComponentStoichiometry(originalComplex, component);
			if (component.equals(componentToReplace)) {
				newComplex.addComponent(replacement);
				if (componentStoichiometry != null) {
					Stoichiometry newComponentStoichiometry = model.addNew(Stoichiometry.class, componentStoichiometry
							.getRDFId().replaceAll(model.getXmlBase(), "") + "_" + splitcount);
					newComponentStoichiometry.setPhysicalEntity(componentToReplace);
					newComponentStoichiometry.setStoichiometricCoefficient(componentStoichiometry
							.getStoichiometricCoefficient());
					splitcount++;
					newElements.add(newComponentStoichiometry);
				}
			} else {
				newComplex.addComponent(component);
				if (componentStoichiometry != null) {
					Stoichiometry newComponentStoichiometry = model.addNew(Stoichiometry.class, componentStoichiometry
							.getRDFId().replaceAll(model.getXmlBase(), "") + "_" + splitcount);
					newComponentStoichiometry.setPhysicalEntity(component);
					newComponentStoichiometry.setStoichiometricCoefficient(componentStoichiometry
							.getStoichiometricCoefficient());
					splitcount++;
					newElements.add(newComponentStoichiometry);
				}
			}
		}
		newComplex.setCellularLocation(originalComplex.getCellularLocation());
		for (Xref xref : originalComplex.getXref()) {
			newComplex.addXref(xref);
		}
		for (String comment : originalComplex.getComment()) {
			newComplex.addComment(comment);
		}
		StringBuilder sb = new StringBuilder();
		for (PhysicalEntity pe : newComplex.getComponent()) {
			sb.append(pe.getDisplayName());
			sb.append(":");
		}
		sb.deleteCharAt(sb.length() - 1);
		newComplex.setDisplayName(sb.toString());
		for (Complex c : model.getObjects(Complex.class)) {
			if (c.equals(newComplex)) {
				continue;
			}
			if (Util.testComplexEquivalency(c, newComplex)) {
				newComplex = c;
				break;
			}
		}
		return newComplex;  //To change body of created methods use File | Settings | File Templates.
	}

	private static boolean testComplexEquivalency(Complex c, Complex newComplex) {
		boolean complexmatch = false;
		if (c.getCellularLocation() != newComplex.getCellularLocation()) {
			if (c.getComponent().equals(newComplex.getComponent()) && c.getComponentStoichiometry().size() ==
					newComplex.getComponentStoichiometry().size()) {
				boolean stoichmatch = true;
				for (Stoichiometry newStoich : newComplex.getComponentStoichiometry()) {
					boolean match = false;
					for (Stoichiometry oldStoich : c.getComponentStoichiometry()) {
						//Float comparison ok as these are manually set/should be consistent - check this pre-launch.
						if (newStoich.getPhysicalEntity().equals(oldStoich.getPhysicalEntity()) &&
								newStoich.getStoichiometricCoefficient() == oldStoich.getStoichiometricCoefficient()) {
							match = true;
							break;
						}
					}
					if (!match) {
						stoichmatch = false;
						break;
					}
				}
				return stoichmatch;
			}
		}
		return complexmatch;
	}

	public static void replaceInControl(Control duplicate, PhysicalEntity toReplace, PhysicalEntity replacement) {
		if (duplicate instanceof Catalysis) {
			Catalysis catalysis = (Catalysis) duplicate;
			HashSet<PhysicalEntity> cofactors = new HashSet<PhysicalEntity>();
			cofactors.addAll(catalysis.getCofactor());
			for (PhysicalEntity cofactor : cofactors) {
				if (cofactor.equals(toReplace)) {
					catalysis.addCofactor(replacement);
					catalysis.removeCofactor(toReplace);
				}
			}
		}
		HashSet<PhysicalEntity> controllerSet = new HashSet<PhysicalEntity>();
		for (Controller ctrlr : duplicate.getController()) {
			if (ctrlr instanceof PhysicalEntity) {
				controllerSet.add((PhysicalEntity) ctrlr);
			}
		}
		for (PhysicalEntity controller : controllerSet) {
			if (controller.equals(toReplace)) {
				duplicate.addController(replacement);
				duplicate.removeController(toReplace);
			}
		}
	}

	public static boolean hasEquivalentInComplex(Complex complex, PhysicalEntity entity) {
		for (PhysicalEntity component : complex.getComponent()) {
			if (component.equals(entity)) {
				return true;
			} else if (component instanceof Protein && entity instanceof Protein) {
				if (component.getComment().contains(KeyParam.REACTOME_ENTITYSET_FLAG) && entity.getComment().contains
						(KeyParam.REACTOME_ENTITYSET_FLAG)) {
					if (component.getMemberPhysicalEntity().size() == entity.getMemberPhysicalEntity().size()) {
						return true;
					}
				} else if (component.getComment().contains(KeyParam.REACTOME_ENTITYSET_FLAG) || entity.getComment()
						.contains(KeyParam.REACTOME_ENTITYSET_FLAG)) {
					Protein set = component.getComment().contains(KeyParam.REACTOME_ENTITYSET_FLAG) ? (Protein)
							component : (Protein) entity;
					Protein single = !component.getComment().contains(KeyParam.REACTOME_ENTITYSET_FLAG) ? (Protein)
							component : (Protein) entity;
					if (set.getMemberPhysicalEntity().contains(single)) {
						return true;
					}
				} else {
					return ((Protein) component).getEntityReference().equals(((Protein) entity).getEntityReference());
				}
			} else {
				return false;
			}
		}
		return false;
	}

	public static PhysicalEntity duplicateWithTransport(Model model, PhysicalEntity physicalEntity, PhysicalEntity
			rightHandSideEntity, long expansionCount) {
		if (physicalEntity instanceof Complex) {
			Complex originalComplex = (Complex) physicalEntity;
			Complex newComplex = model.addNew(Complex.class, rightHandSideEntity.getRDFId().replaceAll(model
					.getXmlBase(), "") + "_" + expansionCount);
			newComplex.setCellularLocation(rightHandSideEntity.getCellularLocation());
			for (PhysicalEntity component : originalComplex.getComponent()) {
				Stoichiometry componentStoichiometry = Util.getComponentStoichiometry(originalComplex, component);
				newComplex.addComponent(duplicateWithTransport(model, component, newComplex, expansionCount));
				expansionCount++;
				if (componentStoichiometry != null) {
					Stoichiometry newComponentStoichiometry = model.addNew(Stoichiometry.class, componentStoichiometry
							.getRDFId().replaceAll(model.getXmlBase(), "") + "_" + expansionCount);
					expansionCount++;
					newComponentStoichiometry.setPhysicalEntity(component);
					newComponentStoichiometry.setStoichiometricCoefficient(componentStoichiometry
							.getStoichiometricCoefficient());
				}
			}
			for (Xref xref : rightHandSideEntity.getXref()) {
				newComplex.addXref(xref);
			}
			for (String comment : rightHandSideEntity.getComment()) {
				newComplex.addComment(comment);
			}
			StringBuilder sb = new StringBuilder();
			for (PhysicalEntity pe : newComplex.getComponent()) {
				sb.append(pe.getDisplayName());
				sb.append(":");
			}
			sb.deleteCharAt(sb.length() - 1);
			newComplex.setDisplayName(sb.toString());
			copyMetadata(rightHandSideEntity, newComplex);
			for (Complex c : model.getObjects(Complex.class)) {
				if (c.equals(newComplex)) {
					continue;
				}
				if (Util.testComplexEquivalency(c, newComplex)) {
					newComplex = c;
					break;
				}
			}
			return newComplex;
		} else if (physicalEntity instanceof Protein) {
			Protein originalProtein = (Protein) physicalEntity;
			Protein newProtein = model.addNew(Protein.class, rightHandSideEntity.getRDFId().replaceAll(model
					.getXmlBase(), "") + "_" + expansionCount);
			expansionCount++;
			newProtein.setCellularLocation(rightHandSideEntity.getCellularLocation());
			for (Xref xref : rightHandSideEntity.getXref()) {
				newProtein.addXref(xref);
			}
			for (String comment : rightHandSideEntity.getComment()) {
				newProtein.addComment(comment);
			}
			newProtein.setDisplayName(rightHandSideEntity.getDisplayName());
			newProtein.setEntityReference(originalProtein.getEntityReference());
			copyMetadata(rightHandSideEntity, newProtein);
			for (Protein p : model.getObjects(Protein.class)) {
				if (p.equals(newProtein)) {
					continue;
				}
				if (p.isEquivalent(newProtein)) {
					newProtein = p;
					break;
				}
			}
			return newProtein;
		} else if (physicalEntity instanceof SmallMolecule) {
			SmallMolecule originalSM = (SmallMolecule) physicalEntity;
			SmallMolecule newSM = model.addNew(SmallMolecule.class, rightHandSideEntity.getRDFId().replaceAll(model
					.getXmlBase(), "") + "_" + expansionCount);
			expansionCount++;
			newSM.setCellularLocation(rightHandSideEntity.getCellularLocation());
			for (Xref xref : rightHandSideEntity.getXref()) {
				newSM.addXref(xref);
			}
			for (String comment : rightHandSideEntity.getComment()) {
				newSM.addComment(comment);
			}
			newSM.setDisplayName(rightHandSideEntity.getDisplayName());
			newSM.setEntityReference(originalSM.getEntityReference());
			copyMetadata(rightHandSideEntity, newSM);
			for (SmallMolecule sm : model.getObjects(SmallMolecule.class)) {
				if (sm.equals(newSM)) {
					continue;
				}
				if (sm.isEquivalent(newSM)) {
					newSM = sm;
					break;
				}
			}
			return newSM;
		} else if (physicalEntity instanceof Rna) {
			Rna originalRNA = (Rna) physicalEntity;
			Rna newRNA = model.addNew(Rna.class, rightHandSideEntity.getRDFId().replaceAll(model.getXmlBase(), "") +
					"_" + expansionCount);
			expansionCount++;
			newRNA.setCellularLocation(rightHandSideEntity.getCellularLocation());
			for (Xref xref : rightHandSideEntity.getXref()) {
				newRNA.addXref(xref);
			}
			for (String comment : rightHandSideEntity.getComment()) {
				newRNA.addComment(comment);
			}
			newRNA.setDisplayName(rightHandSideEntity.getDisplayName());
			newRNA.setEntityReference(originalRNA.getEntityReference());
			copyMetadata(rightHandSideEntity, newRNA);
			for (Rna rna : model.getObjects(Rna.class)) {
				if (rna.equals(newRNA)) {
					continue;
				}
				if (rna.isEquivalent(newRNA)) {
					newRNA = rna;
					break;
				}
			}
			return newRNA;
		} else if (physicalEntity.getClass().equals(PhysicalEntityImpl.class)) {
			PhysicalEntity originalPE = physicalEntity;
			PhysicalEntity newPE = model.addNew(PhysicalEntity.class, rightHandSideEntity.getRDFId().replaceAll(model
					.getXmlBase(), "") + "_" + expansionCount);
			expansionCount++;
			newPE.setCellularLocation(rightHandSideEntity.getCellularLocation());
			for (Xref xref : rightHandSideEntity.getXref()) {
				newPE.addXref(xref);
			}
			for (String comment : rightHandSideEntity.getComment()) {
				newPE.addComment(comment);
			}
			newPE.setDisplayName(rightHandSideEntity.getDisplayName());
			copyMetadata(rightHandSideEntity, newPE);
			for (PhysicalEntity pe : model.getObjects(PhysicalEntity.class)) {
				if (pe.equals(newPE)) {
					continue;
				}
				if (pe.isEquivalent(newPE)) {
					newPE = pe;
					break;
				}
			}
			return newPE;
		} else {
			log.error("Transport didn't catch " + physicalEntity.getRDFId());

		}
		//To change body of created methods use File | Settings | File Templates.
		return null;
	}

	public static ArrayList<Integer> getDifference(double[] control, double[] condition, ProblemInstance
												   instance) {
		ArrayList<Integer> resultList = new ArrayList<Integer>(3);
		for (int i = 0; i < 3; i++) {
			resultList.add(-1);
		}
		try {
			HashSet<Integer> difference = new HashSet<Integer>();
			for (Integer reactionIndex : instance.getReactionSet()) {
				//Test for difference.
				if (control[reactionIndex] != condition[reactionIndex]) {
					difference.add(reactionIndex);
				}
			}
			HashMap<Integer, Object> ioMap  = instance.getIntegerObjectHashMap();
			for (Integer i : difference) {
				Entity e = (Entity)ioMap.get(i);
				System.out.print(e.getRDFId() + "\t");
			}
			System.out.println();
			resultList.set(0, difference.size());
			resultList.set(1, NetworkOutputHandler.generatePredictionOutputSets(difference, instance, true).size());
			resultList.set(2, NetworkOutputHandler.generatePredictionOutputSets(difference, instance, false).size());
		} catch (Exception e) {
			System.out.println("ERROR CONT: " + control.length);
			System.out.println("ERROR COND: " + condition.length);
			e.printStackTrace();
			System.exit(-1);
		}
		return resultList;
	}
}
