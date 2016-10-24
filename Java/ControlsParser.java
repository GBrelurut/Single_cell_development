import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.BufferedWriter;

import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.util.List;
import java.util.ArrayList;

public class ControlsParser {
	public static void main(String[] args) {

		/* Simulate file importation, setting parameters, output and input */
		String annotationFilePath = "/users/csb/brelurut/Documents/Clean_scripts/1_Design_n_Control_Files_preprocessing/mm10ens81.gff";
		String controltableFilePath = "/users/csb/brelurut/Documents/Clean_scripts/1_Design_n_Control_Files_preprocessing/controlstable.tsv";
		File annotationFile = new File(annotationFilePath);
		File controltableFile = new File(controltableFilePath);
		String mitochondrial = "MT";
		String spikes = "ERCC";
		String feature = "gene";
		String attribute = "ID";

		/* Call parsing and writing function */
		writeTable(annotationFile, controltableFile, mitochondrial, spikes,
				feature, attribute);
	}

	/*
	  Write controls table function : 
	  @param : annotation file 
	  @param : output file
	  @param : pattern for mitochondrial genome
	  @param : spike-ins RNA
	  @param : type of feature
	  @param : attribute to consider 
	  @return : void 
	  @purpose : Write a tsv file gathering mitochondrial features and spikes features with their origin
	 */

	public static void writeTable(File annotation, File output, String mtch,
			String spikes, String feature, String attribute) {

		/* Construct patterns from arguments */
		Pattern mt = Pattern.compile(mtch);
		Pattern sp = Pattern.compile(spikes);
		Pattern commentary = Pattern.compile("^#");
		Pattern ft = Pattern.compile(feature);
		Pattern at = Pattern.compile(attribute);

		/* Construct headers of output file */
		String headers = "ID\tType";

		/* Initialise control features lists and line variable*/
		List<String> mitochondrials = new ArrayList<String>();
		List<String> spikeins = new ArrayList<String>();
		String aLine;
		
		/* Parse file and write output */
		try (BufferedReader in = new BufferedReader(new FileReader(annotation));
				BufferedWriter out = new BufferedWriter(new FileWriter(output,
						true))) {

			while ((aLine = in.readLine()) != null) {
				// Do nothing until commentaries end
				Matcher comm = commentary.matcher(aLine);
				if (!comm.find()) {
					// Split fields and set matchers
					String[] parts = aLine.split("\t");
					Matcher mito = mt.matcher(parts[0]);
					Matcher spik = sp.matcher(parts[0]); // this supposes spikes
															// are annotated
															// with first field
															// corresponding to
															// the pattern
					Matcher feat = ft.matcher(parts[2]);
					// Select features else go through
					if (feat.find()) {
						// Extract data for mitochondrial genes
						if (mito.find()) {
							String[] meta = parts[8].split(";");
							for (String data : meta) {
								Matcher attr = at.matcher(data);
								if (attr.find()) {
									String sub = attribute + "=";
									data = data.replaceAll(sub, "");
									mitochondrials.add(data);
									break;
								}
							}
						// Extract data for spikes' genes
						} else if (spik.find()) {
							String[] meta = parts[8].split(";");
							for (String data : meta) {
								Matcher attr = at.matcher(data);
								if (attr.find()) {
									String sub = attribute + "=";
									data = data.replaceAll(sub, "");
									spikeins.add(data);
									break;
								}
							}
						}
					}
				}
			}
			
			// Write results in output file
			out.write(headers);
			for (String el : mitochondrials) {
				out.newLine();
				out.write(el + "\t" + "mitochondrial");
			}
			for (String el : spikeins) {
				out.newLine();
				out.write(el + "\t" + "spike");
			}

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
