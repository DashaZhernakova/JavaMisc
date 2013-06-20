/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.String;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author dashazhernakova
 */
public class SashaJoin {
    public static void main(String[] args) throws IOException {
        BufferedReader sasha_f = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/sasha/Als_COPD_UC_contr_expression_sample_list.txt"));
        BufferedReader lude_f = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/sasha/BloodSATVATLiverMuscleHT12ProbesCentered.txt.25PCAsOverSamplesRemoved.txt.TriTyperFormat.txt"));
        BufferedWriter joint_f = new BufferedWriter(new FileWriter("/Users/dashazhernakova/Documents/UMCG/sasha/Als_COPD_UC_contr_expression_sample_list_res_java.txt"));
        Map <String, String[]> sasha = new HashMap <String, String[]>();
        String[] sasha_header=sasha_f.readLine().split("\t");
        String line = "";
        while ((line=sasha_f.readLine()) != null){
            sasha.put(line.split("\t")[0], line.trim().split("\t"));
        }
        
        String ages = "\t" + "\t" + "\t" + "\t" + "\t" + "\t" + "\t" + "\t";
        String group = "\t" + "\t" + "\t" + "\t" + "\t" + "\t" + "\t" + "\t";
        String[] header = lude_f.readLine().trim().split("\t");
        ArrayList<Integer> leave_pos = new ArrayList<Integer>();
        for (int i = 0; i < header.length; i++){
            String cur_field = header[i];
            if (i < 9)
                joint_f.write(cur_field + "\t");
            else{
                if (sasha.containsKey(cur_field)){
                    leave_pos.add(i);
                    joint_f.write("\t" + cur_field);
                    
                    ages += "\t" + sasha.get(cur_field)[1];
                    group += "\t" + sasha.get(cur_field)[2];
                } 
            }
        }
        joint_f.write("\n" + ages + "\n" + group);    
        int j=0;
        sasha_f.close();
        
        while ((line=lude_f.readLine()) != null){
            joint_f.write("\n");
            String[] line_spl = line.trim().split("\t");
            for (int i = 0; i < 9; i++)
                joint_f.write(line_spl[i] + "\t");
            for (Integer pos : leave_pos)
                joint_f.write("\t" + line_spl[pos]);
            if (j == 10)
                break;
        }
        joint_f.close();
       
    }
    
}
