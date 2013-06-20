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
import java.util.Collections;
import java.util.TreeMap;

/**
 *
 * @author dashazhernakova
 */
public class ProcessSam {
    TreeMap<String, Integer> lengths;
    public ArrayList<String> getFileNames(){
        ArrayList<String> file_names = new ArrayList();
        //File dir = new File("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/");
        File dir = new File("/target/gpfs2/gcc/home/dasha/hapmapData/E-MTAB-197");
        for (File child : dir.listFiles()) {
            if (child.getName().endsWith(".sam")){
                       File f = new File(child.getPath());
                       file_names.add(f.getPath());
            }
        }
        return file_names;
    }
    
    
    
    public ArrayList getAllTranscripts() throws FileNotFoundException, IOException{
        /**
         * get all transcripts from Ensembl
         */
        BufferedReader in = null;
        lengths = new TreeMap();
        ArrayList<String> all_transcr = new ArrayList();
        String line = "";
        //in = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/exmples/summary2.txt"));
        in = new BufferedReader(new FileReader("/target/gpfs2/gcc/home/dasha/GeneticalGenomicsDatasets/hapMapAnnotation_Gnames.txt"));
        //in = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/hapMapAnnotation_Gnames.txt"));
        line = in.readLine();
        while ((line = in.readLine()) != null){
            lengths.put(line.split("\t")[1], (Integer.valueOf(line.split("\t")[5]) - Integer.valueOf(line.split("\t")[4])));
            //lengths.put(line.split("\t")[1], (Integer.valueOf(line.split("\t")[10]) - Integer.valueOf(line.split("\t")[9])+1));
            all_transcr.add(line.split("\t")[1]);
        }
        return all_transcr;
    }
    public void printTable(String[][] table) throws IOException{
        BufferedWriter out_rpkm = new BufferedWriter (new FileWriter("/target/gpfs2/gcc/home/dasha/GeneticalGenomicsDatasets/expression_table_sam_rpkm.txt"));
        BufferedWriter out_reads = new BufferedWriter (new FileWriter("/target/gpfs2/gcc/home/dasha/GeneticalGenomicsDatasets/expression_table_sam_reads.txt"));
        //BufferedWriter out_rpkm = new BufferedWriter (new FileWriter("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/exmples/example_table_sam.txt"));
        //BufferedWriter out_reads = new BufferedWriter (new FileWriter("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/exmples/example_table_sam_reads.txt"));
        String cur_count="", num_of_reads="",tr_name="";
        float rpkm = 0,reads=0;

        for (int tr = 0; tr < table.length; tr++){
            
            for (int sam = 0; sam < table[tr].length; sam++){
                cur_count = table[tr][sam];
                num_of_reads = table[1][sam];
                tr_name=table[tr][0];
                //if (tr_name.equals("ENST00000341421"))
                //        System.out.println("ENST00000341421");
                if (tr == 0)
                    if (sam == 0){
                        out_rpkm.write("");
                        out_reads.write("");
                  
                    }
                    else{ 
                        out_rpkm.write("\t" + cur_count);
                        out_reads.write("\t" + cur_count);
                    }
                else{
                if (sam == 0)
                    //if transcript present in samples
                    if (! cur_count.isEmpty()){

                        out_rpkm.write(cur_count);
                        out_reads.write(cur_count);
                    }
                    else
                        break; //if absent - not write it
                else{
                    rpkm = (Float.valueOf(cur_count)*1000000000)/(Float.valueOf(num_of_reads)*lengths.get(tr_name));
                    //Float.valueOf(table[tr][sam]);
                    reads = (Float.valueOf(cur_count)*1000000)/Float.valueOf(num_of_reads);
                    out_rpkm.write("\t" + rpkm);//RPKM
                    out_reads.write("\t" + reads);//reads
                    //out_reads.write("\t" + reads + "\t" + cur_count + "\t" + num_of_reads+"\t"+lengths.get(tr_name));//reads
                }
                }
            }
            if ((! tr_name.isEmpty()) || (tr == 0)){
                out_rpkm.write("\n");
                out_reads.write("\n");
            }
        }
        out_rpkm.close();
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
    public void readSam() throws IOException{
        BufferedReader in = null;
        
	String line = "";
        String[] spl_line;
        //boolean[] usedTranscr;
        ArrayList<String> all_transcr = new ArrayList();
        //ArrayList acceptedCodes = new ArrayList();
        int index = 0;
	//ArrayList<ArrayList> table = new ArrayList();
        
        
        
        //all transcripts from ensembl, sorted
        all_transcr = getAllTranscripts();
        Collections.sort(all_transcr);
        //usedTranscr = new boolean[all_transcr.size()];
        
        ArrayList<String> fNames = getFileNames();
        //1st row - sample ids; 2nd row - #reads per sample
        //1st column - transcript id
        String[][] table = new String[all_transcr.size()+2][fNames.size() + 1];
        //String[] samples_ids = new String[fNames.length];
        table=fill(table);//fill with empty strings
        
        int cur_id = 1;
        
        for (String f_name : fNames){
            int tr_count=0,r_count=0;
            //ArrayList<String> cur_sample = new ArrayList();
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
               
                if (! line.startsWith("@")){
                    spl_line = line.split("\t");
                    String qual = spl_line[4];
                    String rname = spl_line[2];
                    String mate_rname = spl_line[6]; 
                    //get the transcript position in all transcripts list 
                    index = Collections.binarySearch(all_transcr, rname);
      
                    if ((! rname.equals("*")) && (Integer.valueOf(qual)>=10) && (mate_rname.equals("="))){

                          //check
                          
                              if (index >= 0){
                          index = index + 2;//to leave the first row in the table for sample ids
                              r_count++;
                              //if there are reads mapped to the transcript, add 1 to the #reads(tr, sample) 
                              if (!table[index][cur_id].equals("0.0")){
                                table[index][cur_id] = Integer.toString(Integer.valueOf(table[index][cur_id].trim()).intValue() + 1);       
                                table[1][cur_id] = Integer.toString(Integer.valueOf(table[1][cur_id].trim()).intValue() + 1);      
                              }
                              //if a new transcript 
                              else{
                                tr_count++;
                                table[index][0] = rname;//write transcript id in the forst column. if empty, this transcript doesn't occur in the mapping                                  
                                if (table[0][cur_id].equals("")){
                                    table[0][cur_id]=sample_id; //write sample id in the first row
                                    table[1][cur_id]="1"; //write #reads per sample
                                    table[index][cur_id] = "1"; //write # reads
                                }
                                else{
                                    table[index][cur_id] = "1"; //write # reads
                                    table[1][cur_id] = Integer.toString(Integer.valueOf(table[1][cur_id].trim()).intValue() + 1);
                                }
                                   
                             }
                        }
                    }
                }
            }
            cur_id++;
            System.out.println("tr_count="+tr_count+" ;r_count=" + r_count);
        }
        //print the result
        printTable(table);
    }
    
    public static void main(String[] args) throws IOException {

        ProcessSam p = new ProcessSam();
        p.readSam();
        
    }
}
