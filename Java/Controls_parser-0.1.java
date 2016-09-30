import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.BufferedWriter;

import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.util.List;
import java.util.ArrayList;

public class Controls_parser{
	public static void main(String [] args){
		/* Simulate file importation, setting parameters, output and input */
		String annotationFilePath="/users/csb/brelurut/Documents/Clean_scripts/1_Design_n_Control_Files_preprocessing/mm10ens81.gff";
		String controltableFilePath= "/users/csb/brelurut/Documents/Clean_scripts/1_Design_n_Control_Files_preprocessing/controlstable.tsv";
		File annotationFile=new File (annotationFilePath);
		File controltableFile=new File (controltableFilePath);
		String mitochondrial="MT";
		String spikes = "ERCC";
		String feature ="gene";
		String attribute="ID";
		
		/* Call parsing and writing function */
		write_table(annotationFile, controltableFile, mitochondrial, spikes, feature, attribute);
	}
	
	/* Write controls table function :
	 input : annotation file
	 others arguments : pattern for mitochondrial genome, spike-ins RNA, type of feature, attribute to consider
	 output : void
	 purpose : Write a tsv file gathering mitochondrial features and spikes features with their origin */
	
	public static void write_table(File annotation, File output, String mtch, String spikes, String feature, String attribute){
		
		
		/* Construct patterns from arguments */
		Pattern mt = Pattern.compile(mtch);
		Pattern sp = Pattern.compile(spikes);
		Pattern commentary = Pattern.compile("^#");
		Pattern ft=Pattern.compile(feature);
		Pattern at= Pattern.compile(attribute);
		
		/* Construct headers of output file */
		String headers = "ID\tType";
		
		/* Initialize writer */
		FileWriter fstream = null;
		BufferedWriter out = null;
		
		/* Initialize control features lists */
		List <String> mitochondrials = new ArrayList<String>();
		List <String> spikeins = new ArrayList<String>();
		
		/* Initialize input for reading */
		FileInputStream fis;
		
		/* Parse file seeking on corresponding fields :
			a) feature pattern
			b) mitochondrial or spikes pattern
			c) attribute pattern on split attributes field
			then modifying attribute removing "*=" */
		try {
			fis = new FileInputStream(annotation);
			BufferedReader in = new BufferedReader(new InputStreamReader(fis));
			String aLine;
			while ((aLine = in.readLine()) != null){
				Matcher comm = commentary.matcher(aLine);
				if(! comm.find()){
					String [] parts= aLine.split("\t");
					Matcher mito = mt.matcher(parts[0]);
					Matcher spik = sp.matcher(parts[0]); // this suppose spikes are annotated with first field corresponding to the pattern
					Matcher feat = ft.matcher(parts[2]);
					if (feat.find()){
						if(mito.find()){
							String[] meta = parts[8].split(";"); 
							for(String data : meta){
								Matcher attr = at.matcher(data);
								if(attr.find()){
									String sub= attribute +"=";
									data = data.replaceAll(sub, "") ;
									mitochondrials.add(data);
									break;
								}
							}
						}
						else if(spik.find()){
							String[] meta = parts[8].split(";"); 
							for(String data : meta){
								Matcher attr = at.matcher(data);
								if(attr.find()){
									String sub= attribute +"=";
									data = data.replaceAll(sub, "") ;
									spikeins.add(data);
									break;
								}
							}
						}
					}
				}
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
		/* Write output file */
		try {
			fstream = new FileWriter(output, true);
			out = new BufferedWriter(fstream);
		} catch(IOException e1){
			e1.printStackTrace();
		}
			
		try{
			out.write(headers);
			for(String el : mitochondrials){
				out.newLine();
				out.write(el +"\t" + "mitochondrial");
			}
			for(String el : spikeins){
				out.newLine();
				out.write(el + "\t" + "spike");
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		try {
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}