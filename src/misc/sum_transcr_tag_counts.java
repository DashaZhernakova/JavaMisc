/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author dashazhernakova
 */
public class sum_transcr_tag_counts {
    float[] tot; //total number of tags per individual
    String header;
    TreeMap<String, List<String>> enst2lines;
    TreeSet<String> tags2excl;
    
    /**
     * old method
     * @throws IOException 
     */
    public void changeTagNameToTranscrNames() throws IOException{
        BufferedReader covbed_out = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/ForLude_tagwise_hg19/res_out.txt"));
        BufferedReader expr_file = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/ForLude_tagwise_hg19/tagwise_expression_table.txt"));
        //BufferedReader covbed_out = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/tmp_cov.txt"));
        //BufferedReader expr_file = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/tmp.txt"));
        
        BufferedWriter out = new BufferedWriter(new FileWriter("/Users/dashazhernakova/Documents/UMCG/ForLude_tagwise_hg19/exonwise_expression.txt"));
        
        
        String[] spl;
        String line="";
        Map<String, String[]> tags = new HashMap<String, String[]>();
        while ((line = covbed_out.readLine()) != null){
            spl = line.split("\t");
            tags.put(spl[0], spl);
        }
        line=expr_file.readLine();   
        out.write(line + "\n");
        
        String[] cur_tag;
        while ((line = expr_file.readLine()) != null){
            spl = line.split("\t");
            cur_tag = tags.get(spl[0].replaceAll("_.?$", "").replace("chr", ""));
            if (cur_tag != null){
                //System.out.println(line);
                //System.out.println(cur_tag);
                //out.write(cur_tag[1] + "\t" + cur_tag[0]);
                out.write(cur_tag[1].replaceFirst(":.*", ""));
                for (int i = 1; i < spl.length; i++)
                    out.write("\t" + Float.toString(Float.valueOf(spl[i]) * Float.valueOf(cur_tag[2])));
                out.write("\n");
            }
        }
        
       out.close();
        
    }
    
    /**
     * fills a map enst2lines: transcript : list of all tagwise lines for this transcript
     * @throws IOException 
     */
    public void changeTagNameToTranscrNamesSimple(String tagwise_expression, String tag2transcripts) throws IOException{
        BufferedReader covbed_out = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/ForLude_tagwise_hg19/res_out_uniq.txt"));
        //BufferedReader expr_file = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/deepSAGE_tag/tagwise_expression_table_SNP_in_recognition_sequence_tags_excluded.txt"));
        TextFile tags2transcr_file = new TextFile(tag2transcripts, false);
        TextFile expr_file = new TextFile(tagwise_expression, false);
        expr_file.open();
        //BufferedReader covbed_out = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/tmp_cov.txt"));
        //BufferedReader expr_file = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/tmp.txt"));

        String[] spl;
        String line="";
        
        Map<String, String[]> tags = new HashMap<String, String[]>();
        Map<String, String> tag2transcr = new HashMap<String, String>();
        /*while ((line = covbed_out.readLine()) != null){
            spl = line.split("\t");
            tags.put(spl[0], spl);
        }
         * 
         */
        while ((spl = tags2transcr_file.readLineElems(TextFile.tab)) != null)
            tag2transcr.put(spl[0], spl[1]);
        
        line=expr_file.readLine();   
        header = line;
        
        int len = line.split("\t").length;
        float[] totalPerSample = new float[len];
        
        enst2lines = new TreeMap<String, List<String>>();
        //String[] cur_tag;
        String cur_tr;
        while ((line = expr_file.readLine()) != null){
            spl = line.split("\t");
            //cur_tag = tags.get(spl[0].replaceAll("_.?$", "").replace("chr", ""));
            cur_tr = tag2transcr.get(spl[0]);
            if (cur_tr != null){
                String cur_line = "";
                for (int i = 1; i < len; i++){
                    cur_line += String.format("\t%.2f", Float.valueOf(spl[i]));
                    totalPerSample[i] += Float.valueOf(spl[i]);
                }
                
                //adding a line with tag-wise expression to the transcript entry in enst2lines
                //String cur_tr = cur_tag[1].replaceFirst(":.*", "");
                //for exon level summarization:
                
                //String cur_tr = cur_tag[1];
                
                List<String> linesThisFar = enst2lines.get(cur_tr);
                if (linesThisFar == null) {
                    linesThisFar = new ArrayList<String>();
                    linesThisFar.add(cur_line);
                    enst2lines.put(cur_tr, linesThisFar);
                }
                else
                    enst2lines.get(cur_tr).add(cur_line);

            }
        }
        //total numeber of tags per individual
        tot = totalPerSample;
    }
    
    /**
     * sum up expression values for transcripts in enst2lines
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public void sumUpSimpleMap(String output) throws FileNotFoundException, IOException{
        //BufferedWriter out = new BufferedWriter(new FileWriter("/Users/dashazhernakova/Documents/UMCG/ForLude_tagwise_hg19/exonwise_expression_table_SNP_in_recognition_sequence_tags_excluded.txt"));
        TextFile out = new TextFile(output, true);
        List<String> lines = new ArrayList<String>();
        //BufferedReader in = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/tmp_res.txt"));
  
        out.write(header);
        int l = header.split("\t").length;
        for (String tr : enst2lines.keySet()){
            lines = enst2lines.get(tr);
            float [] cur_sum = new float[l];
            for (String line : lines){
                String[] cur_spl = line.split("\t");
                for (int i = 1; i < l; i++)
                    cur_sum[i]+=Float.valueOf(cur_spl[i]);
            }
            out.write("\n" + tr);
            for (int i = 1; i < l; i++){        
                out.write("\t" + String.format("%.5f", cur_sum[i]));
                //out.write("\t" + String.format("%.4f", cur_sum[i]*1000000/tot[i]));
            }       
        }
        out.close();
    }
    
    /**
     * excludes tags from tags2exc
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public void excludeTags() throws FileNotFoundException, IOException{
        BufferedReader tags2exc = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/ForLude_tagwise_hg19/Tags_with_SNP_in_recognition_sequence.txt"));
        String line = "";
        tags2excl = new TreeSet<String>();
        TreeSet<String> tmp = new TreeSet<String>();
        while ( (line = tags2exc.readLine()) != null ){
            tags2excl.add(line.split("\t")[0].replace("chr", "").replaceAll("_[/+/-]$", ""));
        }
        
        BufferedReader tag_expr = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/deepSAGE_tag/tagwise_expression_table.txt"));
        BufferedWriter out = new BufferedWriter(new FileWriter("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/deepSAGE_tag/tagwise_expression_table_SNP_in_recognition_sequence_tags_excluded.txt"));
        
        while ( (line = tag_expr.readLine()) != null ){
             if ( ! tags2excl.contains(line.split("\t")[0]))
                out.write(line + "\n");
        }        
        out.close();      
    }
    
    /**
     * old method
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public void sumUp() throws FileNotFoundException, IOException{
        BufferedReader in = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/deepSAGE_tag/tagwise_expression_table_sorted.txt"));
        //BufferedReader in = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/deepSAGE_tag/tmp.txt"));
        BufferedWriter out = new BufferedWriter(new FileWriter("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/deepSAGE_tag/tagwise_expression_table_sorted_summed.txt"));
        
        
        out.write(in.readLine());
        String s = in.readLine();
        String[] cur_spl = s.split("\t");
        String cur_tr = cur_spl[0];
        int l = cur_spl.length - 1;
        float [] cur_sum = new float[l];
        while (s!=null){
            for(int i=1;i<=l;i++)
                cur_sum[i-1]+=Float.valueOf(cur_spl[i]);            
            s = in.readLine();
            if (s==null){
                cur_spl =new String[1]; cur_spl[0]="Null";
            }else 
                cur_spl = s.split("\t");
            if (!cur_spl[0].equals(cur_tr)){
                //sum and writeln
                out.write("\n"+cur_tr);
                for(int i=0;i<l;i++)
                    out.write("\t"+Float.toString(cur_sum[i]));
                //init again
                cur_tr = cur_spl[0];
                cur_sum = new float[l];
            }
        }
        out.close();
    
    }
    
    /**
     * old method
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public void sumUpSimple() throws FileNotFoundException, IOException{
    BufferedReader in = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/ForLude_tagwise_hg19/transcriptwise_expression_sorted_simple.txt"));
    BufferedWriter out = new BufferedWriter(new FileWriter("/Users/dashazhernakova/Documents/UMCG/ForLude_tagwise_hg19/res_transcriptwise_expression_simple.txt"));
        
    //BufferedReader in = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/tmp_res.txt"));
    //BufferedWriter out = new BufferedWriter(new FileWriter("/Users/dashazhernakova/Documents/tmp_res_su.txt"));
    out.write(in.readLine());
        //out.write(header);
        String s = in.readLine();
        String[] cur_spl = s.split("\t");
        String cur_tr = cur_spl[0];
        int l = cur_spl.length - 1;
        float [] cur_sum = new float[l];
        while (s!=null){
            for(int i=1;i<=l;i++)
                cur_sum[i-1]+=Float.valueOf(cur_spl[i]);            
            s = in.readLine();
            if (s==null){
                cur_spl =new String[1]; cur_spl[0]="Null";
            }else 
                cur_spl = s.split("\t");
            if (!cur_spl[0].equals(cur_tr)){
                //sum and writeln
                out.write("\n"+cur_tr);
                for(int i=0;i<l;i++){
                    out.write("\t"+Double.toString(cur_sum[i]));
                    //System.out.println(cur_tr + " : " + tot_loc[i+1]);
                    
                }
                //init again
                cur_tr = cur_spl[0];
                cur_sum = new float[l];
            }
        }
        out.close();
    
    }
    
    public static void main(String[] args) throws FileNotFoundException, IOException{
        /*BufferedReader b = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/chip-seq/expression_transcriptwise_sorted.txt"));
        BufferedWriter w = new BufferedWriter(new FileWriter("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/chip-seq/expression_transcriptwise_res.txt"));
        w.write(b.readLine());
        String s = b.readLine();
        String[] cur_spl = s.split("\t");
        String ct = cur_spl[0];
        int l = cur_spl.length - 1;
        int [] csum = new int[l];
        while (s!=null){
            for(int i=1;i<=l;i++)csum[i-1]+=Integer.valueOf(cur_spl[i]);            
            s = b.readLine();
            if (s==null){
                cur_spl =new String[1]; cur_spl[0]="Null";
            }else cur_spl = s.split("\t");
            if (!cur_spl[0].equals(ct)){
                //sum and writeln
                w.write("\n"+ct);
                for(int i=0;i<l;i++)w.write("\t"+Integer.toString(csum[i]));
                //init again
                ct = cur_spl[0];
                csum = new int[l];
            }
        }
        w.close();
    }*/
        sum_transcr_tag_counts s = new sum_transcr_tag_counts();
        //s.excludeTags();
        s.changeTagNameToTranscrNamesSimple("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/deepSAGE_tag/tagwise_expression_table_SNP_in_recognition_sequence_tags_excluded_NORMALIZED.txt",
                "/Users/dashazhernakova/Documents/UMCG/ForLude_tagwise_hg19/tag_positions_NoSNPsInRecognitionSeq_to_lincRNAs.txt");
        s.sumUpSimpleMap("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/deepSAGE_transcr/NORMALIZED/lincRNA/expresion_table.txt");
        //s.sumUp();
        //s.excludeTags();
        //System.out.println("asdf_fsfdg_-".replaceAll("_[/+/-]$", ""));
        
        
    }
}
