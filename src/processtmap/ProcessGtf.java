/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package processtmap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

/**
 *
 * @author Dasha
 */
public class ProcessGtf {

    /**
     * get an expressionData table for eQTL mapping pipeline from cuffcompare .tmap file
     */
    public ArrayList getAllTranscripts() throws FileNotFoundException, IOException{
        /**
         * get all transcripts from Ensembl
         */
        BufferedReader in = null;
        ArrayList all_transcr = new ArrayList();
        String line = "";
        in = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/hapMapAnnotation_Gnames.txt"));
        //in = new BufferedReader(new FileReader("/target/gpfs2/gcc/home/dasha/GeneticalGenomicsDatasets/hapMapAnnotation_Gnames.txt"));
        while ((line = in.readLine()) != null)
            all_transcr.add(line.split("\t")[1]);
        return all_transcr;
    }
    public ArrayList<String> getFileNames(){
        ArrayList<String> file_names = new ArrayList();
        //File dir = new File("/target/gpfs2/gcc/home/dasha/tophat_out/emtab1");
        /*File dir = new File("/Users/dashazhernakova/Documents/UMCG/cluster/");
        for (File ch : dir.listFiles()) {
            if (ch.isDirectory())
                for (File child : ch.listFiles()){
                   if (child.getName().endsWith("cufflinks_out_bwa")){
                       File f = new File(child.getPath() +"/transcripts.gtf");
                       if (f.exists())
                        file_names.add(f.getPath());
                  }
                }
        }*/
        file_names.add("/Users/dashazhernakova/Documents/UMCG/hapmapData/tophat/1234_5.gtf");
        return file_names;
    }
    
    public void printTable(String[][] table) throws IOException{
        //BufferedWriter out = new BufferedWriter (new FileWriter("/target/gpfs2/gcc/home/dasha/GeneticalGenomicsDatasets/bwaToCuff_expression_table.txt"));
        BufferedWriter out = new BufferedWriter (new FileWriter("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/bwa+cuff.txt"));
        String cur_count="", num_of_reads="",tr_name="";

        for (int tr = 0; tr < table.length; tr++){
            
            for (int sam = 0; sam < table[tr].length; sam++){
                cur_count = table[tr][sam];
                num_of_reads = table[1][sam];
                tr_name=table[tr][0];
                if (tr == 0)
                    if (sam == 0) out.write("");
                    else out.write("\t" + cur_count);
                else{
                if (sam == 0)
                    //if transcript present in samples
                    if (! cur_count.isEmpty())
                        out.write(cur_count);
                    else
                        break; //if absent - not write it
                else
                    out.write("\t" + cur_count);//RPKM
                }
            }
            if ((! tr_name.isEmpty()) || (tr == 0))
                out.write("\n");
        }
        out.close();
    }
    
    public String[][] fill(String[][] table){
        for (int i = 0; i < table.length; i++)
            for (int j = 0; j < table[i].length; j++)
                if ((i == 0) || (j ==0))
                    table[i][j] = "";
                else
                    table[i][j] = "0.0";
        return table;
    }
    
    public void readTmap() throws IOException{
        BufferedReader in = null;       
	String line = "";
        String[] spl_line;
        ArrayList<String> all_transcr = new ArrayList();
        int index = 0;
        boolean skip=false;       
        //all transcripts from ensembl, sorted
        all_transcr = getAllTranscripts();
        Collections.sort(all_transcr);
 
        ArrayList<String> fNames = getFileNames();
        String[][] table = new String[all_transcr.size()+1][fNames.size() + 1];

        table=fill(table);//fill with empty strings
        String type="", tr_name="",cur_count="";
        float rpkm=0;
        int cur_id = 1;
        
        for (String f_name : fNames){
            in = new BufferedReader(new FileReader(f_name));
            //getting sample name from the file name
            System.out.println("processing file: " + f_name);
            String sample_id="";
            for (String sampleid : f_name.split("/")){
            if (! sampleid.isEmpty())
               if (Character.isDigit(sampleid.charAt(0))){
                   if (sampleid.length() > 10) skip = true;
                   sample_id = sampleid.substring(0, 6);
                   break;
               }
            }
            if (skip)        
                continue;
            while ((line = in.readLine()) != null){
               spl_line = line.split("\t");
               type=spl_line[2]; //transcript or exon
               tr_name = spl_line[0];

               if (type.equals("transcript")){                             
                    //get the transcript position in all transcripts list 
                    index = Collections.binarySearch(all_transcr, tr_name);
                    //if found in all transcripts
                    if (index >= 0){
                        rpkm = Float.valueOf(spl_line[8].split(";")[2].split("\"")[1]).floatValue();
                        //read_count = Float.valueOf(spl_line[7].trim()).floatValue()*Float.valueOf(spl_line[12].trim()).floatValue()/1000; 
                        index = index + 1;//to leave the first row in the table for sample ids
                        cur_count = table[index][cur_id];        
                          //if there are reads mapped to the transcript, summ up  RPKMs 
                         if (!cur_count.equals("0.0"))
                            table[index][cur_id] = Float.toString(Float.valueOf(cur_count.trim()).floatValue() + rpkm);      
                         //if a new transcript 
                         else{
                            table[index][0] = tr_name;//write transcript id in the forst column. if empty, this transcript doesn't occur in the mapping                                  
                            table[0][cur_id]=sample_id; //write sample id in the first row
                            table[index][cur_id] = Float.toString(rpkm); //write RPKM
                         }
                    }
                }
            }
            cur_id++;
        }
        //print the result
        printTable(table);
    }
    public static void main(String[] args) throws IOException {

        ProcessGtf p = new ProcessGtf();
        p.readTmap();
        
    }
}
