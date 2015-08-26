package FastDAS;

import IO.BioPAXInputHandler;
import IO.BioPAXOutputHandler;
import IO.TextOutputHandler;
import LogicalSystem.ProblemInstance;
import ModelCuration.BucketFixer;
import ModelCuration.ComplexFlattener;
import ModelCuration.LoopBreaker;
import org.apache.log4j.Level;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.Control;
import org.biopax.paxtools.model.level3.Conversion;
import org.biopax.paxtools.model.level3.PhysicalEntity;

import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: uqlfearn
 * Date: 15/08/13
 * Time: 4:05 PM
 * To change this template use File | Settings | File Templates.
 */
public class MainClass {

	private static org.apache.log4j.Logger log = Logger
			.getLogger(MainClass.class);

	public static void main(String[] argsv) {

		//autoCuratorStart();
	}


	public static void dispatch(File biopaxFile, File pathwayFile, File speciesFile, int objectiveBoxSelectedIndex,
								int optModeComboBoxSelectedIndex, int outputTypeComboBoxSelectedIndex,
								File outputFile, boolean saveMidConf, boolean saveLowConf, String
										experimentDescription) {
		Model model = BioPAXInputHandler.getModel(biopaxFile, pathwayFile, speciesFile);
		System.out.println("INTO WRITE");
		ProblemInstance problemInstance = new ProblemInstance(model, objectiveBoxSelectedIndex,
				optModeComboBoxSelectedIndex, outputTypeComboBoxSelectedIndex, outputFile, saveMidConf, saveLowConf,
				experimentDescription);
		//IO.TextOutputHandler.writeSaveFileForEditing(model, outputFile);
	}

	public static void autoCuratorStart(File owlFile, File pathwayFile, File speciesFile, File outputDirectory) {
		LogManager.getRootLogger().setLevel(Level.ERROR);
		//Read the model into memory
		Model model = BioPAXInputHandler.getModel(owlFile, pathwayFile, speciesFile);
		System.out.println("Model contains : " + (model.getObjects(Conversion.class).size()) + " conversions, "
				+ (model.getObjects(Control.class).size()) + " control elements and " + model.getObjects
				(PhysicalEntity.class).size() + " physical entities");
		long participantcount = 0;
		for (PhysicalEntity pe : model.getObjects(PhysicalEntity.class)) {
			if (!pe.getParticipantOf().isEmpty()) {
				participantcount++;
			}
		}
		System.out.println("Of which " + participantcount + " physicalEntities participate in reactions");
		//Flatten entitySets and Complexes where possible (note that some EntitySets/Complexes will not be flattened
		// by this method).
		ComplexFlattener.flattenEntitySetsInModelFile(model);
		Model m = ComplexFlattener.flattenComplexesInModelFile(model);
		//Write the initial network representation to a SIF representation for visualisation/analysis in Cytoscape.
		TextOutputHandler.writeSIFFromNetwork(m, outputDirectory.getAbsolutePath() + System.getProperty("file" +
				".separator") + "initialNetwork.sif");
		//Perform model debucketing.
		BucketFixer.debucketModel(m, outputDirectory);
		//Write the intermediate (post-debucketing, pre-loopbreaking) network representation to a SIF representation
		// for visualisation/analysis in Cytoscape.
		TextOutputHandler.writeSIFFromNetwork(m, outputDirectory.getAbsolutePath() + System.getProperty("file" +
				".separator") + "debucketNetwork.sif");
		//Perform automatic loop breaking.
		LoopBreaker.breakLoops(m, outputDirectory);
		//Write the curated network representation to a SIF representation for visualisation/analysis in Cytoscape.
		TextOutputHandler.writeSIFFromNetwork(m, outputDirectory.getAbsolutePath() + System.getProperty("file" +
				".separator") + "postCurationNetwork.sif");
		//Output a BioPAX level 3 representation of the model.
		BioPAXOutputHandler.writeBioPaxModel(m, outputDirectory.getAbsolutePath() + System.getProperty("file" +
				".separator") + "curatedModel.owl");
		//Clean the model file (fix any rdf:about in place of rdf:ID, add '#' characters where necessary)
		BioPAXOutputHandler.cleanModel(outputDirectory.getAbsolutePath() + System.getProperty("file.separator"),
				"curatedModel");
		System.out.println("Model contains : " + (m.getObjects(Conversion.class).size()) + " conversions, "
				+ (m.getObjects(Control.class).size()) + " control elements and " + m.getObjects(PhysicalEntity.class)
				.size() + " physical entities");
		participantcount = 0;
		for (PhysicalEntity pe : m.getObjects(PhysicalEntity.class)) {
			if (!pe.getParticipantOf().isEmpty()) {
				participantcount++;
			}
		}
		System.out.println("Of which " + participantcount + " physicalEntities participate in reactions");
	}
}
