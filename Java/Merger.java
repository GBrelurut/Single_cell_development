import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.BufferedWriter;

import java.util.List;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class Merger {
	public static void main(String[] args) throws NullPointerException,
			IOException {
		/* Simulate file importation, setting output and input */
		String srcFile1Path = "/users/csb/brelurut/Documents/Clean_scripts/1_Design_n_Control_Files_preprocessing/readcountstep_output_expression_MDP10.tsv";
		String srcFile2Path = "/users/csb/brelurut/Documents/Clean_scripts/1_Design_n_Control_Files_preprocessing/readcountstep_output_expression_MDP11.tsv";
		String srcFile3Path = "/users/csb/brelurut/Documents/Clean_scripts/1_Design_n_Control_Files_preprocessing/readcountstep_output_expression_MDP12.tsv";
		String mergedFilePath = "/users/csb/brelurut/Documents/Clean_scripts/1_Design_n_Control_Files_preprocessing/mergedFile.tsv";

		/* Construct files and files List */
		List<File> files = new ArrayList<File>();
		files.add(new File(srcFile1Path));
		files.add(new File(srcFile2Path));
		files.add(new File(srcFile3Path));
		File mergedFile = new File(mergedFilePath);

		/* Call merging function */
		mergeFiles(files, mergedFile);
	}

	/*
	  Write controls table function : 
	  @param : list of expression files
	  @param : output file
	  @return : void 
	  @purpose : Write a tsv file gathering expression data
	 */
	
	public static void mergeFiles(List<File> files, File mergedFile)
			throws NullPointerException, IOException {

		/* Test parameters */
		if (files == null) {
			throw new NullPointerException("files argument cannot be null");
		}

		if (mergedFile == null) {
			throw new NullPointerException("mergedFile argument cannot be null");
		}

		/* Set structures and extract first file */
		File init = files.get(0);
		List<String> genes = new ArrayList<>();
		Map<String, String> lines = new HashMap<>();

		// TO DO : write method for genes name extraction
		String names = "MDP10\tMDP11\tMDP12";

		/* Treat first file initialising structures */
		try (BufferedReader in = new BufferedReader(new FileReader(init))) {
			String aLine;
			String[] parts = new String[2];
			int i = 0;
			while ((aLine = in.readLine()) != null) {
				if (i > 0) {
					parts = aLine.split("\t");
					genes.add(parts[0]);
					lines.put(parts[0], parts[1]);
				}
				i++;
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		/* Remove file treated from files list */
		files.remove(0);

		/* Treat other files */
		for (File file : files) {
			try (BufferedReader in = new BufferedReader(new FileReader(file))) {
				String aLine;
				String[] parts = new String[2];
				while ((aLine = in.readLine()) != null) {
					parts = aLine.split("\t");
					String value = lines.get(parts[0]);
					value = value + "\t" + parts[1];
					lines.put(parts[0], value);
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		/* Write merged file, excluding non detected genes */
		try (BufferedWriter out = new BufferedWriter(new FileWriter(mergedFile,
				true))) {
			out.write(names);
			out.newLine();
			for (String key : genes) {
				int sum = 0;
				String[] values = lines.get(key).split("\t");
				for (String value : values) {
					sum = sum + Integer.parseInt(value);
				}
				if (sum > 0) {
					out.write(key + "\t" + lines.get(key));
					out.newLine();
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}