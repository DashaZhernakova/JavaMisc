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
import java.util.HashMap;
import java.util.Map;
/**
 *
 * @author dashazhernakova
 */
public class Compare {
    public void compare() throws FileNotFoundException, IOException{
        BufferedReader bwacuff = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/bwa+cuff.txt"));
        BufferedReader bwa = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromBWA/expression_rpkm.txt"));
        BufferedReader cuff = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromCuff/expression_table_fromCuff.txt"));
        BufferedReader htseq = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/hapmapData/htseq-count_strand.txt"));
        
        String line="", bwacuff_str="",bwa_str="",cuff_str="", tmp_bwa="",tmp_cuff, tmp_htseq="",htseq_str="";
        String[] tr;
        ArrayList<String> all_tr = new ArrayList();
        Map<String,String> bwa_tr = new HashMap<String,String>();
        Map<String,String> cuff_tr = new HashMap<String,String>();
        Map<String,String> htseq_tr = new HashMap<String,String>();
        while ((line = bwa.readLine()) != null){
            if (line.split("\t")[22].startsWith("1382"))
                    System.out.println("1382");
            else
                bwa_tr.put(line.split("\t")[0], line.split("\t")[22]);
        }
        while ((line = cuff.readLine()) != null){
            if (line.split("\t")[1].startsWith("1382"))
                    System.out.println("1382");
            else 
                cuff_tr.put(line.split("\t")[0], line.split("\t")[1]);
        }
        while ((line = htseq.readLine()) != null){
            if (line.split("\t")[1].startsWith("1382"))
                    System.out.println("1382");
            else 
                htseq_tr.put(line.split("\t")[0], line.split("\t")[1]);
        }
        while ((line = bwacuff.readLine()) != null){
            tr=line.split("\t");
            tmp_bwa=bwa_tr.get(tr[0]);
            tmp_cuff =cuff_tr.get(tr[0]);
            tmp_htseq=htseq_tr.get(tr[0]);
            if ((tmp_cuff != null) && (tmp_bwa != null) && (tmp_htseq != null)){
                bwacuff_str+=","+tr[1];
                bwa_str+=","+tmp_bwa;
                cuff_str+=","+tmp_cuff;
                htseq_str+=","+tmp_htseq;
            }
        }
        
        bwacuff_str=bwacuff_str.replaceFirst(",", "bwacuff<-c(")+")";
        bwa_str=bwa_str.replaceFirst(",", "bwa<-c(")+")";
        cuff_str=cuff_str.replaceFirst(",", "cuff<-c(")+")";
        htseq_str=htseq_str.replaceFirst(",", "htseq<-c(")+")";
        System.out.println(htseq_str);
        System.out.println(bwacuff_str);
        System.out.println(bwa_str);
        System.out.println(cuff_str);
    }
    public static void main(String[] args) throws IOException {
        Compare c = new Compare();
        c.compare();
    }
    
}
