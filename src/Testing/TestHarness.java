package Testing;

import FastDAS.KeyParam;
import IO.BioPAXInputHandler;
import LogicalSystem.AdjList;
import LogicalSystem.ProblemInstance;
import LogicalSystem.Util;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.Interaction;
import org.biopax.paxtools.model.level3.PhysicalEntity;
import org.biopax.paxtools.model.level3.Complex;
import org.biopax.paxtools.model.level3.Protein;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Workspace for developing and testing specific methods before moving them into the main project.
 *
 * Created by lfearnley on 5/02/15.
 */
public class TestHarness {

	public static void main (String[] argsv) {
		try {
			Model m = BioPAXInputHandler.readModelFromFile(new File("/Users/lfearnley/Dropbox/LoFParse/LoFparse/curatedModel.owl"));
			HashMap<String, HashSet<String>> uniprotToRDF = new HashMap<String, HashSet<String>>();
			for (Protein p : m.getObjects(Protein.class)) {
				String up = Util.getUniProt(p);
				if (!up.isEmpty()) {
					if (!uniprotToRDF.containsKey(up)) {
						uniprotToRDF.put(up, new HashSet<String>());
					}
					uniprotToRDF.get(up).add(p.getRDFId());
				}
			}
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File("outfile.txt")));
			for (String key : uniprotToRDF.keySet()) {
				bw.write(key);
				bw.write("\t");
				for (String s : uniprotToRDF.get(key)) {
					bw.write(s);
					bw.write(",");
				}
				bw.newLine();
			}
			bw.flush();
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
