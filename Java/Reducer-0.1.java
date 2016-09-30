import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.BufferedWriter;

import java.util.regex.Pattern;
import java.util.regex.Matcher;

public class Reducer {
	public static void main(String[] args){
		/* Simulate file importation, setting output and input */ 
		String designFilePath= "/users/csb/brelurut/Documents/Clean_scripts/1_Design_n_Control_Files_preprocessing/design.txt";
		String phenoFilePath= "/users/csb/brelurut/Documents/Clean_scripts/1_Design_n_Control_Files_preprocessing/phenoData.tsv";
		
		File designFile= new File (designFilePath);
		File phenoFile = new File (phenoFilePath);
		
		design_reduce(designFile, phenoFile);
		
	}
	
	
	/* design reduce function :
	input : design file for Eoulsan
	output : void
	purpose : reads design file and write table in the Columns Section into a phenoData file. */
	
	public static void design_reduce(File file, File output){
		/* Initialize writer */
		FileWriter fstream = null;
		BufferedWriter out = null;
		
		/* Initialize boolean and Columns pattern*/ 
		boolean table = false;
		String pattern = "(Columns)";
		Pattern seek = Pattern.compile(pattern);
		
		/* Initialize input */
		FileInputStream fis;
		
		/* Read file, testing pattern, and writing after pattern is found */
		try {
			fis = new FileInputStream(file);
			BufferedReader in = new BufferedReader(new InputStreamReader(fis));
			String aLine;
			while ((aLine = in.readLine()) != null){
				if(table==true){
					System.out.println(aLine);
					try {
						fstream = new FileWriter(output, true);
						out = new BufferedWriter(fstream);
					} catch(IOException e1){
						e1.printStackTrace();
					}
					try{
						out.write(aLine);
						out.newLine();
					} catch (IOException e) {
						e.printStackTrace();
					}
					try {
						out.close();
					} catch (IOException e) {
						e.printStackTrace();
					}
				}else {
					Matcher m = seek.matcher(aLine);
					if (m.find()){
						table=true;
					}
				}
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
}