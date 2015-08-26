package IO;

import FastDAS.KeyParam;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.Model;

import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created with IntelliJ IDEA.
 * User: uqlfearn
 * Date: 16/08/13
 * Time: 1:04 PM
 * To change this template use File | Settings | File Templates.
 */
public class BioPAXOutputHandler {

	private static Pattern resourcePattern = Pattern.compile("(resource\\=\")[^#]");

	/**
	 * Fairly straightforward method for writing a BioPAX Level 3 format OWL file to disk. Note that this may have
	 * errors in it - see the comments on IO.cleanModel().
	 *
	 * @param m    a paxtools Model object
	 * @param name a String file name for the output model, used in the FileOutputStream constructor.
	 */
	public static void writeBioPaxModel(Model m, String name) {
		try {
			SimpleIOHandler handler;
			handler = new SimpleIOHandler();
			//Uncomment and comment SimpleReader for Jena IO.
			// JenaIOHandler IOHandler = new JenaIOHandler();
			handler.convertToOWL(m, new FileOutputStream(name));
		} catch (Exception e) {
			System.out.println("Exception thrown in BioPax input handling class.");
			e.printStackTrace();
		}
	}

	/**
	 * Current versions of the paxtools library create new objects in a model with rdf:about tags instead of rdf:ID.
	 * This method is invoked to clean this up, replacing all occurrences of rdf:about with rdf:ID. This is not done in
	 * a sophisticated way - it is potentially possible to create a namespace conflict with two identical rdf:IDs in a
	 * file. This also adds in missing pound (#) characters to resource references.
	 *
	 * @param filepath a String containing the path to the model directory
	 * @param filename a String containing the filename for the model being cleaned
	 */
	public static void cleanModel(String filepath, String filename) {
		try {
			File originalModelFile = new File(filepath + filename + ".owl");
			File fixedModelFile = new File(filepath + filename + "_clean.owl");
			BufferedReader br = new BufferedReader(new FileReader(originalModelFile));
			BufferedWriter bw = new BufferedWriter(new FileWriter(fixedModelFile));
			String current = null;
			while ((current = br.readLine()) != null) {
				current = current.replaceAll("rdf:about", "rdf:ID");
				Matcher matcher = resourcePattern.matcher(current);
				if (matcher.find()) {
					current = current.replace(matcher.group(1), "resource=\"#");
				}
				bw.write(current);
				bw.write(" \t  \t\t   \t ");
				bw.write(KeyParam.NEWLINE);
			}
			bw.flush();
			bw.close();
			br.close();
			originalModelFile.delete();
			fixedModelFile.renameTo(originalModelFile);
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

}
