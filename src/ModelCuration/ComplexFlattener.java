package ModelCuration;

import FastDAS.KeyParam;
import LogicalSystem.Util;
import org.apache.log4j.Logger;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.CellularLocationVocabulary;
import org.biopax.paxtools.model.level3.Complex;
import org.biopax.paxtools.model.level3.PhysicalEntity;
import org.biopax.paxtools.model.level3.Stoichiometry;

import java.util.HashSet;
import java.util.concurrent.LinkedBlockingQueue;

/**
 * Created with IntelliJ IDEA.
 * User: uqlfearn
 * Date: 26/08/13
 * Time: 10:32 AM
 * To change this template use File | Settings | File Templates.
 */
public class ComplexFlattener {

	private static org.apache.log4j.Logger log = Logger.getLogger(ComplexFlattener.class);


	public static void flattenEntitySetsInModelFile(Model model) {
		try {
			HashSet<String> affected = new HashSet<String>();
			HashSet<PhysicalEntity> entitySets = new HashSet<PhysicalEntity>();
			for (PhysicalEntity pe : model.getObjects(PhysicalEntity.class)) {
				if (pe.getComment().contains(KeyParam.REACTOME_ENTITYSET_FLAG)) {
					entitySets.add(pe);
				}
			}
			for (PhysicalEntity pe : entitySets) {
				LinkedBlockingQueue<PhysicalEntity> memberQueue = new LinkedBlockingQueue<PhysicalEntity>();
				memberQueue.add(pe);
				while (!memberQueue.isEmpty()) {
					PhysicalEntity currentPE = memberQueue.poll();
					if (currentPE.getComment().contains(KeyParam.REACTOME_ENTITYSET_FLAG)) {
						memberQueue.addAll(currentPE.getMemberPhysicalEntity());
						pe.removeMemberPhysicalEntity(currentPE);
					} else {
						affected.add(pe.getRDFId());
						pe.addMemberPhysicalEntity(currentPE);
					}
				}
			}

		} catch (Exception e) {

		}
	}

	//File modelFile
	public static Model flattenComplexesInModelFile(Model model) {
		//Model model = BioPAXInputHandler.readModelFromFile(modelFile);
		try {
			long stoichCount = 0;
			HashSet<Complex> complexSet = new HashSet<Complex>(model.getObjects(Complex.class));
			//Create a queue to hold modified complexes for collision tests after flattening.
			LinkedBlockingQueue<Complex> modifiedComplexProcessingQueue = new LinkedBlockingQueue<Complex>();
			for (Complex complex : complexSet) {
				//If a complex doesn't have any components, we test to see whether it's an EntitySet. If the Complex
				//is an EntitySet then we continue, otherwise we log it for reporting to Reactome as an error -
				//complexes should always have components.
				if (complex.getComponent() == null || complex.getComponent().isEmpty()) {
					if (complex.getComment().contains(KeyParam.REACTOME_ENTITYSET_FLAG)) {
						continue;
					} else {
						log.warn("WARN : " + complex.getRDFId().replaceAll(model.getXmlBase(), "") + " has no " +
								"components");
						continue;
					}
				}
				//Determine whether the complex needs modification. We use a 'flag' to avoid
				//ConcurrentModificationExceptions (the getComponent method returns the single representation of a
				//Complex's components, not a copy of it).
				boolean process = false;
				for (PhysicalEntity component : complex.getComponent()) {
					if (component instanceof Complex) {
						//Process the complex.
						process = true;
						break;
					}
				}
				if (process) {
					//Add the modified complex to the processing queue for collision resolution.
					modifiedComplexProcessingQueue.add(processComplex(complex, model, stoichCount));
				}
			}
			fixCollisions(modifiedComplexProcessingQueue, model);
			for (Complex c : model.getObjects(Complex.class)) {
				LinkedBlockingQueue<PhysicalEntity> peQueue = new LinkedBlockingQueue<PhysicalEntity>();
				peQueue.addAll((c.getComponent()));
				while (!peQueue.isEmpty()) {
					PhysicalEntity pe = peQueue.poll();
					if (pe.getComment().contains(KeyParam.REACTOME_ENTITYSET_FLAG)) {
						if (BucketFixer.needsExpansion(pe)) {
							c.addComment(KeyParam.REACTOME_COMPLEX_ENTITYSET_CONTAINING_FLAG);
						}
					} else if (pe instanceof Complex) {
						log.error("ERROR: non flat complex after flattening");
						System.exit(-1);
					}
				}
			}
		} catch (Exception e) {
			//model = BioPAXInputHandler.readModelFromFile(modelFile);
		} finally {
			return model;
		}
	}

	/**
	 * We need to iterate through the modified complexes to determine whether any of them are duplicates/collide -
	 * duplication
	 * of entities breaks network connectivity and has caused major issues in models in the past. This method is
	 * written
	 * in such a way that it assumes that other PhysicalEntities within the network are unique, i.e given two
	 * representations
	 * of A (A1, A2), the complexes C1 = A1:B, C2 = A2:B will be treated as separate entities. If we had two complexes
	 * C1 = A1:B, C3 = A1:B and they were in the same physical location, these would be merged.
	 *
	 * @param modifiedComplexProcessingQueue A queue containing Complex objects for testing and processing
	 * @param model                          The model containing the Complex objects
	 */
	private static void fixCollisions(LinkedBlockingQueue<Complex> modifiedComplexProcessingQueue, Model model) {
		try {
			while (!modifiedComplexProcessingQueue.isEmpty()) {
				Complex currentComplex = modifiedComplexProcessingQueue.poll();
				CellularLocationVocabulary currentComplexLocation = currentComplex.getCellularLocation();
				for (Complex c : model.getObjects(Complex.class)) {
					//Test physical locations:
					if (c.getCellularLocation().equals(currentComplexLocation)) {
						//If the Complexes are co-located, test component identities:
						if (c.getComponent() != null && c.getComponent().equals(currentComplex.getComponent())) {
							//Test stoichiometries:
							HashSet<Stoichiometry> stoichSet = new HashSet<Stoichiometry>(currentComplex
									.getComponentStoichiometry());
							for (Stoichiometry stoichiometry : c.getComponentStoichiometry()) {
								stoichSet.remove(Util.getEquivalentStoichiometryFromSet(stoichiometry, stoichSet));
							}
							//If the stoichSet is empty, we've found an equivalent Complex - same location, same
							//components, same stoich of components:
							if (stoichSet.isEmpty()) {
								Util.replaceEntity(currentComplex, c, model);
							}
						}
					}
				}
			}
		} catch (Exception e) {

		} finally {

		}
	}

	/**
	 * Flattens an individual complex, returning the modified complex. This method is invoked when a Complex contains
	 * one or more Complex objects amongst its components. If this happens, we need to 'flatten' the Complex, i.e.,
	 * shift
	 * components up from subComplexes. This is an iterative process - complexes should contain only non-complexes as
	 * components. This is one exception to this - if a Complex object contains EntitySet Complexes, then we have to
	 * leave them in-place.
	 *
	 * @param complex A Complex containing other Complexes as its components.
	 * @param model   The Model containing the Complex
	 * @return the flattened Complex
	 */
	private static Complex processComplex(Complex complex, Model model, long stoichCount) {
		try {
			//Create a queue to hold unprocessed components
			LinkedBlockingQueue<PhysicalEntity> componentQueue = new LinkedBlockingQueue<PhysicalEntity>();
			//Initialise queue to contain all components of the complex
			componentQueue.addAll(complex.getComponent());
			while (!componentQueue.isEmpty()) {
				PhysicalEntity current = componentQueue.poll();
				//If the current component is an instance of a Complex, and the component complex contains other
				// components
				//(i.e. is not an EntitySet or an error in curation), we want to replace the component complex with its
				//components (recursively, so we're adding non-Complex PhysicalEntity objects, not more Complexes):
				if (current instanceof Complex && !((Complex) current).getComponent().isEmpty()) {
					Complex subComplex = (Complex) current;
					//Process the subComplex. This is necessary because the subComplex may have Complexes in its
					// components
					processComplex(subComplex, model, stoichCount);
					//Get the Stoichiometry associated with the subComplex. This will return null if there's no
					// stoichiometry
					//associated with the subComplex.
					Stoichiometry correspondingStoich = Util.getComponentStoichiometry(complex, subComplex);
					float subComplexOriginalComplex = (float) 1.0;
					if (correspondingStoich != null) {
						subComplexOriginalComplex = correspondingStoich.getStoichiometricCoefficient();
					}
					//For each component of the subComplex:
					for (PhysicalEntity pe : subComplex.getComponent()) {
						//Sort out stoichiometry.
						//Let s(a, b) be the stoichiometry of some entity a in the complex b. Then, the stoichiometry
						//of some entity x in the final fixedComplex given some originalComplex and its subComplex is:
						//s(x, fixedComplex) = s(subComplex, originalComplex)*s(x, subComplex)+s(x, originalComplex)
						float xSubComplex = (float) 1.0;
						Stoichiometry xSubComplexStoich = Util.getComponentStoichiometry(subComplex, pe);
						if (xSubComplexStoich != null) {
							xSubComplex = xSubComplexStoich.getStoichiometricCoefficient();
						}
						//Assume that there's no pe in the parent complex - this makes xOriginalComplex = 0.
						float xOriginalComplex = (float) 0.0;
						//And if there's no pe the stoichiometry object should pass back a null.
						Stoichiometry xOriginalComplexStoich = null;
						//If the complex doesn't have the pe already, add it.
						if (!complex.getComponent().contains(pe)) {
							complex.addComponent(pe);
							//If the complex has the pe, then the xOriginalComplex is either 1.0 (if no stoich
							// specified) or
							//the specified stoich value:
						} else {
							xOriginalComplexStoich = Util.getComponentStoichiometry(complex, pe);
							if (xOriginalComplexStoich != null) {
								xOriginalComplex = xOriginalComplexStoich.getStoichiometricCoefficient();
							} else {
								xOriginalComplex = (float) 1.0;
							}
						}
						float xFixedComplex = subComplexOriginalComplex * xSubComplex + xOriginalComplex;
						//If we're updating a Stoich object, do so now:
						if (xOriginalComplexStoich != null) {
							xOriginalComplexStoich.setStoichiometricCoefficient(xFixedComplex);
							//Otherwise we need to make one if the xFixedComplex value != 1.0:
						} else if (Math.round(xFixedComplex) != 1) {
							Stoichiometry newStoich = model.addNew(Stoichiometry.class,
									correspondingStoich.getRDFId().replaceAll(model.getXmlBase(), "") + "_EXPAND" +
											stoichCount);
							stoichCount++;
							newStoich.setStoichiometricCoefficient(xFixedComplex);
							newStoich.setPhysicalEntity(pe);
							complex.addComponentStoichiometry(newStoich);
						}
					}
					//We now remove the subComplex stoichiometry (if any) from the Complex
					if (correspondingStoich != null) {
						complex.removeComponentStoichiometry(correspondingStoich);
						//We can also remove the stoichiometry from the Model - the component that it referred to will
						//be removed from the Complex, so it will otherwise be an orphan.
						model.remove(correspondingStoich);
					}
					//Remove the subComplex from the Complex
					complex.removeComponent(subComplex);
					//We don't need to alter anything that isn't a Complex - this gets left in place.
				} else {
					continue;
				}
			}
		} catch (Exception e) {
			log.error("ERROR : Something has failed in the Complex flattener: ", e);
			//TODO : Remove or comment out System.exit before upload/publication.
			System.exit(-1);
		} finally {
			return complex;
		}
	}
}
