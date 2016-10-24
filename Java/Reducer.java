import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.BufferedWriter;

import java.util.regex.Pattern;
import java.util.regex.Matcher;

public class Reducer {
	public static void main(String[] args) throws IOException {
		/* Simulate file importation, setting output and input */
		String designFilePath = "/users/csb/brelurut/Documents/Clean_scripts/1_Design_n_Control_Files_preprocessing/design.txt";
		String phenoFilePath = "/users/csb/brelurut/Documents/Clean_scripts/1_Design_n_Control_Files_preprocessing/phenoData.tsv";
		File designFile = new File(designFilePath);
		File phenoFile = new File(phenoFilePath);

		/* Call reducing function */
		designReduce(designFile, phenoFile);

	}

	/*
	  Write controls table function : 
	  @param : design file
	  @param : output file
	  @return : void 
	  @purpose : Write a tsv file gathering samples metadata
	 */
	
	public static void designReduce(File file, File output) throws IOException {
		/* Initialise variables */
		String aLine;
		boolean table = false;
		String pattern = "(Columns)";
		Pattern seek = Pattern.compile(pattern);

		/* Read file, copy lines */
		try(BufferedReader in = new BufferedReader(new FileReader(file))) {		
			while ((aLine = in.readLine()) != null) {
				if (table == true) {
					try(BufferedWriter out = new BufferedWriter(new FileWriter(output, true)) ) {
						out.write(aLine);
						out.newLine();
					} catch (IOException e2) {
						e2.printStackTrace();
					}
				} else {
					Matcher m = seek.matcher(aLine);
					if (m.find()) {
						table = true;
					}
				}
			}
		} catch (IOException e1) {
			e1.printStackTrace();
		}

	}
}