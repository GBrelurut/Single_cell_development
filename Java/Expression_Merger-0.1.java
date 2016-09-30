import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.BufferedWriter;

import java.util.List;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class Merger {
	public static void main(String[] args){
		/* Simulate file importation, setting output and input */
		String srcFile1Path = "/users/csb/brelurut/Documents/Clean_scripts/1_Design_n_Control_Files_preprocessing/readcountstep_output_expression_MDP10.tsv";
		String srcFile2Path = "/users/csb/brelurut/Documents/Clean_scripts/1_Design_n_Control_Files_preprocessing/readcountstep_output_expression_MDP11.tsv";
		String srcFile3Path = "/users/csb/brelurut/Documents/Clean_scripts/1_Design_n_Control_Files_preprocessing/readcountstep_output_expression_MDP12.tsv";
		
		String mergedFilePath= "/users/csb/brelurut/Documents/Clean_scripts/1_Design_n_Control_Files_preprocessing/mergedFile.tsv";
		
		/* Construct files and files List */
		List<File> files = new ArrayList<File>();
		files.add(new File(srcFile1Path));
		files.add(new File(srcFile2Path));
		files.add(new File(srcFile3Path));
		File mergedFile = new File(mergedFilePath);
		
		/* Call Merging Function */
		mergeFiles(files, mergedFile);
	}
	
	/* Merge Files function :
	 input : List of Files
	 output : void
	 purpose : Write a file merging lines of the input files, first line is the sample names */
	 
	public static void mergeFiles(List <File> files, File mergedFile){
		/* Initialize writer */
		FileWriter fstream = null;
		BufferedWriter out = null;
		
		/* Initialize variable container */
		List<String> genes = new ArrayList<String>();
		Map<String,String> lines = new HashMap<String,String>();
		
		/* Simulate first line construction */
		String names = "MDP10\tMDP11\tMDP12";
		
		/* First iteration, initiate values of List and Map */
		File init = files.get(0);
		FileInputStream fis;
		try {
			fis = new FileInputStream(init);
			BufferedReader in = new BufferedReader(new InputStreamReader(fis));
			String aLine;
			String[] parts = new String[2];
			int i = 0;
			while ((aLine = in.readLine()) != null){
				if(i>0){
					parts= aLine.split("\t");
					genes.add(parts[0]);
					lines.put(parts[0], parts[1]);
				}
				i++;
			}
			in.close();
		} catch (IOException e) { // IOException is not dispensable
			e.printStackTrace();
		}
		/* Remove first file to iterate through following files */
		files.remove(0);
		for( File file : files){
			try {
				fis = new FileInputStream(file);
				BufferedReader in = new BufferedReader(new InputStreamReader(fis));
				String aLine;
				int i = 0;
				String[] parts = new String[2];
				while ((aLine = in.readLine()) != null){
					parts= aLine.split("\t");
					String value = lines.get(parts[0]);
					value = value + "\t" + parts[1];
					lines.put(parts[0],value);
				}
				in.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		/* Create output file */
		try {
			fstream = new FileWriter(mergedFile, true);
			out = new BufferedWriter(fstream);
		} catch(IOException e1){
			e1.printStackTrace();
		}
			
		/* Write first line */
		try{
			out.write(names);
			out.newLine();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		/* Write following lines */
		for(String key : genes){
			int sum = 0;
			String [] values = lines.get(key).split("\t");
			for(String value : values){
				sum = sum + Integer.parseInt(value);
			}
			if(sum > 0){
				try{
					out.write(key+"\t"+lines.get(key));
					out.newLine();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		try {
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}