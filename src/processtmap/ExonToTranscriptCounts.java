package processtmap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

public class ExonToTranscriptCounts {
    
    public void process(String f_name) throws FileNotFoundException, IOException{
        BufferedReader ex_counts=new BufferedReader(new FileReader(f_name));
        BufferedWriter tr_counts=new BufferedWriter(new FileWriter(f_name.replace(".txt", "_out.txt")));
        String line="", cur_transcr="", count="", transcr="";
        int cur_count=0;
        while ((line = ex_counts.readLine()) != null){
            String[] fields = line.split("\t");
            count = fields[0];
            transcr = fields[1];
            if (transcr.equals(cur_transcr))
                cur_count += Integer.valueOf(count);
            else{
                if (!cur_transcr.isEmpty()) 
                    tr_counts.write(cur_transcr + "\t" + Integer.toString(cur_count)+"\n");
                cur_count=Integer.valueOf(count);
                cur_transcr=transcr;
            }
        }
        tr_counts.write(cur_transcr + "\t" + Integer.toString(cur_count));
                
        tr_counts.close();
        ex_counts.close();       
    }
    public static void main(String[] args) throws IOException {
        System.out.println(args);
        ExonToTranscriptCounts p = new ExonToTranscriptCounts();
        p.process("/Users/dashazhernakova/Documents/UMCG/cluster/tmp.txt");
        //p.process("/target/gpfs2/gcc/home/dasha/tophat_out/emtab1/1382_1_1/tmp.txt");
    }
}
