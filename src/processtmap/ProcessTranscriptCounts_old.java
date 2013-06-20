/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package processtmap;

import java.io.BufferedReader;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

/**
 * processes files of the type:
 * tr_id\tcount
 * tr_id\tcount
 * 
 * generates a table of ids vs tr_counts
 * ...
 * @author dashazhernakova
 */
public class ProcessTranscriptCounts_old {
    
    String annotation_file_name;
    String dir_name;
    
    ProcessTranscriptCounts_old(String dir, String annot){
        annotation_file_name=annot;
        dir_name=dir;
        
    }
    public ArrayList getAllTranscripts() throws FileNotFoundException, IOException{
        /**
         * get all transcripts from Ensembl
         */
        BufferedReader in = null;
        ArrayList all_transcr = new ArrayList();
        String line = "";
        in = new BufferedReader(new FileReader(annotation_file_name));
        //in = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/hapMapAnnotation_Gnames.txt"));
        //in = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/cluster/scripts/test_all_tr.txt"));
        //in = new BufferedReader(new FileReader("/target/gpfs2/gcc/home/dasha/GeneticalGenomicsDatasets/hapMapAnnotation_Gnames.txt"));
        while ((line = in.readLine()) != null)
            all_transcr.add(line.split("\t")[1]);
        return all_transcr;
    }
    public ArrayList<String> getFileNames(){
        ArrayList<String> file_names = new ArrayList();
        //File dir = new File("/target/gpfs2/gcc/home/dasha/tophat_out/emtab1");
        //File dir = new File("/Users/dashazhernakova/Documents/UMCG/hapmapData/tophat/");
        File dir = new File(dir_name);
        for (File ch : dir.listFiles()) {
            if (ch.isDirectory())
                for (File child : ch.listFiles()){
                   if (child.getName().endsWith("exon_counts_sorted_out.txt")){
                        file_names.add(child.getPath());
                  }
            }
        }
        //file_names.add("/Users/dashazhernakova/Documents/UMCG/cluster/scripts/test1.txt");
        //file_names.add("/Users/dashazhernakova/Documents/UMCG/cluster/scripts/test2.txt");
        //file_names.add("/Users/dashazhernakova/Documents/UMCG/hapmapData/tophat/1671_5_1/htseq_out.txt");
        return file_names;
    }
    
    public void printTable(String[][] table) throws IOException{
        //BufferedWriter out_reads = new BufferedWriter (new FileWriter("/target/gpfs2/gcc/home/dasha/GeneticalGenomicsDatasets/expression_table_sam_reads.txt"));
        BufferedWriter out_reads = new BufferedWriter (new FileWriter("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_htseq/res_table_stranded.txt"));
        String cur_count="", num_of_reads="",tr_name="";
        float reads=0;

        for (int tr = 0; tr < table.length; tr++){
            
            for (int sam = 0; sam < table[tr].length; sam++){
                cur_count = table[tr][sam];
                num_of_reads = table[1][sam];
                tr_name=table[tr][0];
                //if (tr_name.equals("ENST00000341421"))
                //        System.out.println("ENST00000341421");
                if (tr == 0)
                    if (sam == 0){
                        out_reads.write("");
                  
                    }
                    else{ 
                        out_reads.write("\t" + cur_count);
                    }
                else{
                if (sam == 0)
                    //if transcript present in samples
                    if (! cur_count.isEmpty()){
                        out_reads.write(cur_count);
                    }
                    else
                        break; //if absent - not write it
                else{
                    reads = (Float.valueOf(cur_count)*1000000)/Float.valueOf(num_of_reads);
                    out_reads.write("\t" + reads);//reads
                    //out_reads.write("\t" + reads + "\t" + cur_count + "\t" + num_of_reads+"\t"+lengths.get(tr_name));//reads
                }
                }
            }
            if ((! tr_name.isEmpty()) || (tr == 0)){
                out_reads.write("\n");
            }
        }
        out_reads.close();
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
    public void readCounts() throws IOException{
        BufferedReader in = null;
        
	String line = "";
        String[] spl_line;
        ArrayList<String> all_transcr = new ArrayList();
        int index = 0;
	
        
        //all transcripts from ensembl, sorted
        all_transcr = getAllTranscripts();
        Collections.sort(all_transcr);
        
        ArrayList<String> fNames = getFileNames();
        //1st row - sample ids; 2nd row - #reads per sample
        //1st column - transcript id
        String[][] table = new String[all_transcr.size()+2][fNames.size() + 1];
        
        table=fill(table);//fill with empty strings
        
        int cur_id = 1;
        
        for (String f_name : fNames){
            int tr_count=0,r_count=0;
            in = new BufferedReader(new FileReader(f_name));
            //getting sample name from the file name
            System.out.println("processing file: " + f_name);
            String sample_id="";
            for (String sampleid : f_name.split("/")){
            if (! sampleid.isEmpty())
               if (Character.isDigit(sampleid.charAt(0))){
                   sample_id = sampleid.substring(0, 6);
                   break;
               }
            }
            
                    
            while ((line = in.readLine()) != null){
                    spl_line = line.split("\t");
                    String tr_id = spl_line[0];
                    String count = spl_line[1];
                    
                    //get the transcript position in all transcripts list 
                    index = Collections.binarySearch(all_transcr, tr_id);
                           
                    if (index >= 0){ //if transcript found
                        index = index + 2;//to leave the first row in the table for sample ids
                        r_count+=Integer.valueOf(count).intValue();
                        
                        //if there are reads mapped to the transcript, add count to the #reads(tr, sample) 
                        if (!table[index][cur_id].equals("0.0")){
                            //sum up #reads for this sample and this transcript
                            table[index][cur_id] = Integer.toString(Integer.valueOf(table[index][cur_id].trim()).intValue() + Integer.valueOf(count).intValue());       
                            //add count to the overall number of reads for this sample
                            table[1][cur_id] = Integer.toString(Integer.valueOf(table[1][cur_id].trim()).intValue() + Integer.valueOf(count).intValue());      
                        }
                        //if a new transcript 
                        else{
                            tr_count++;
                            table[index][0] = tr_id;//write transcript id in the first column. if empty, this transcript doesn't occur in the mapping                                  
                            if (table[0][cur_id].equals("")){//if that's the first entry of the sample 
                                table[0][cur_id]=sample_id; //write sample id in the first row
                                table[1][cur_id]=count; //write #reads per sample
                                table[index][cur_id] = count; //write # reads
                            }
                            else{
                                table[index][cur_id] = count; //write # reads
                                table[1][cur_id] = Integer.toString(Integer.valueOf(table[1][cur_id].trim()).intValue() + Integer.valueOf(count).intValue());//add #reads per sample
                            }
                                   
                        }
                   }
                    
                
            }
            cur_id++;
            System.out.println("overall number of reads=" + r_count);
        }
        //print the result
        printTable(table);
    }
    public static void main(String[] args) throws IOException {

        ProcessTranscriptCounts_old p;
        if (args.length > 1){
            p = new ProcessTranscriptCounts_old(args[0], args[1]);
            p.readCounts();
        }
        else
            System.out.println("no input file given");
        
    }
}


