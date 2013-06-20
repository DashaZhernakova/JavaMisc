package misc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 * combines 2 expression tables: if value present in both - average it, if only in one - write it.
 * If individuals don't correspond to each other exits. Manually insert the absent individuals with expression values filled with '*'
 * @author dashazhernakova
 */

public class CombineYaleArgonne {
    HashSet<String> yUnique; // transcripts that appear only in Yale
    HashSet<String> aUnique; // transcripts that appear only in Argonne
    ArrayList<String> yMixups; // list of mixups in Yale
    ArrayList<String> aMixups; // list of mixups in Argonne
    ArrayList<Integer> yMixupsPos; // list of mixups in Yale
    ArrayList<Integer> aMixupsPos; // list of mixups in Argonne
    HashSet<Integer> y_dupl = new HashSet<Integer>();
    HashSet<Integer> a_dupl = new HashSet<Integer>();
    int numLines;
    public CombineYaleArgonne(){
     yUnique = new HashSet<String>();
     aUnique = new HashSet<String>();
     
     yMixups = new ArrayList<String>();
     aMixups = new ArrayList<String>();
     yMixupsPos = new ArrayList<Integer>();
     aMixupsPos = new ArrayList<Integer>();
     
     y_dupl = new HashSet<Integer>();
     a_dupl = new HashSet<Integer>();
        
     numLines = 0;
}
    public CombineYaleArgonne(String aMx, String yMx){
     yUnique = new HashSet<String>();
     aUnique = new HashSet<String>();
     
     yMixups = new ArrayList<String>();
     aMixups = new ArrayList<String>();
     yMixupsPos = new ArrayList<Integer>();
     aMixupsPos = new ArrayList<Integer>();
     
     yMixups.addAll(Arrays.asList(yMx.split(";")));
     aMixups.addAll(Arrays.asList(aMx.split(";")));
     
     y_dupl = new HashSet<Integer>();
     a_dupl = new HashSet<Integer>();
     
     numLines = 0;
     
}
    /**
     * checks if individuals are the same for the 2 files
     * @param header1
     * @param header2
     * @return true if individuals are the same, false - otherwise
     */
    public boolean checkSamplePositions(String header1, String header2){       
        String[] spl1 = header1.replaceAll("_yale", "").split("\t");
        String[] spl2 = header2.replaceAll("_argonne", "").split("\t");

        ArrayList<String> h1 = new ArrayList<String>();
        h1.addAll(Arrays.asList(spl1));
        
        for (int i = 1; i < spl1.length; i ++){
            if (yMixups.contains(spl1[i]))
                yMixupsPos.add(i);
            if (spl1[i].matches("NA[0-9]*_[0-9].*"))
                y_dupl.add(i);

            
        }
        for (int i = 1; i < spl2.length; i ++){
            if (aMixups.contains(spl2[i]))
                aMixupsPos.add(i);
            if (spl2[i].matches("NA[0-9]*_[0-9].*"))
                y_dupl.add(i);
        }
            //System.out.println(spl1[i] + " " + spl2[i]);    
            //if (! spl1[i].equals(spl2[i]))
              //      System.out.println(spl1[i] + " " + spl2[i]);
        if (spl1.length != spl2.length){ 
            System.out.println("Length of the header differs");
            return false;
        }
        else
            for (int i = 1; i < spl1.length; i ++)
                if (! spl1[i].equals(spl2[i])){
                    System.out.println(spl1[i] + " != " + spl2[i]);
                    return false;
                }
        
        return true;
    }
    
    
    
    /**
     * finds unique transcripts for both files
     * @param yale_f - file name 1
     * @param arg_f - file name 2
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public void findUnique(String yale_f, String arg_f) throws FileNotFoundException, IOException{
        BufferedReader yale = new BufferedReader(new FileReader(yale_f));
        BufferedReader arg = new BufferedReader(new FileReader(arg_f));
        TextFile yu = new TextFile("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Yale_non-coding/y_unique", true);
        TextFile au = new TextFile("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Yale_non-coding/a_unique", true);
        yu.open();
        au.open();
        
        String y = yale.readLine();
        String a = arg.readLine();
        HashSet<String> yList = new HashSet<String>();
        HashSet<String> aList = new HashSet<String>();
        int a_lines = 0;
        while ((y = yale.readLine()) != null){
             yList.add(y.split("\t")[0]);
             numLines++;
         }
        while ((a = arg.readLine()) != null){
             aList.add(a.split("\t")[0]);
             a_lines++;
         }
        if (a_lines > numLines)
            numLines = a_lines;
        for (String tr : yList)
            if (!aList.contains(tr)){
                yu.writeln(tr);
                yUnique.add(tr);
           }
        for (String tr : aList)
            if (!yList.contains(tr)){
                au.writeln(tr);
                aUnique.add(tr);
            }
         System.out.println("unique probes " + aUnique.size()+yUnique.size());
         yale.close();
         arg.close();
         
         yu.close();
         au.close();
         
    }
    /**
     * combines 2 tables. If individuals don't correspond to each other exits. 
     * ! Manually insert the absent individuals with expression values filled with '*'
     * @param yale_f - file name
     * @param arg_f
     * @param f_name - output
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public void combine(String yale_f, String arg_f, String f_name) throws FileNotFoundException, IOException{
        BufferedReader yale = new BufferedReader(new FileReader(yale_f));
        BufferedReader arg = new BufferedReader(new FileReader(arg_f));
        BufferedWriter out = new BufferedWriter(new FileWriter(f_name));
        
        System.out.println("looking for unique probes....");
        findUnique(yale_f, arg_f);
       
        String y = yale.readLine();
        String a = arg.readLine();
        String comb_line = "";
        int cur_num = 0;
        Float cur_avg;
        int l = a.split("\t").length;
        System.out.println("Done\nlooking for sample mess....");
        if (! checkSamplePositions(y,a)){
            System.out.println("columns are not identical!!!");
            return;
        }
        System.out.println("started combining tables");
        out.write(a.replaceAll("_argonne", ""));
        
        y = yale.readLine();
        a = arg.readLine();
        while ((y != null) && (a != null)){
            if (cur_num % 10000 == 0)
                System.out.println("processed " + cur_num + " of " + numLines);
            cur_num++; 
            String[] y_spl = y.split("\t");
            String[] a_spl = a.split("\t");
            
            //if transcript ids differ
            if (!a_spl[0].equals(y_spl[0])){    
                if (aUnique.contains(a_spl[0])){
                    out.write("\n" + a);
                    a = arg.readLine();
                    
                }
                if (yUnique.contains(y_spl[0])){
                    out.write("\n" + y);
                    y = yale.readLine();
                }
            }
            //if the same transcript ids
            else{
                comb_line = a_spl[0];
                for (int ind = 1 ; ind < l; ind++){
                    String cur_y = y_spl[ind];
                    String cur_a = a_spl[ind];
                    
                    if ( ( cur_y.equals("*") && (cur_a.equals("*") ) ||
                            (yMixupsPos.contains(ind)) && (aMixupsPos.contains(ind)) ) ) 
                        System.out.println("NB! a mixup in both datasets! : " + cur_a + "," + cur_y);
                    
                    else{ //if not a mixup in both
                        if ( (cur_y.equals("*")) || (yMixupsPos.contains(ind)) ) //to skip mixups and absent samples
                            cur_avg = Float.valueOf(cur_a);
                        else if ( (cur_a.equals("*")) || (aMixupsPos.contains(ind)) ) //to skip mixups and absent samples
                            cur_avg = Float.valueOf(cur_y);
                        else
                            cur_avg = (Float.valueOf(cur_y) + Float.valueOf(cur_a))/2; //averaging expression
                        comb_line += "\t" + Float.toString(cur_avg);
                    }
                }
            
                out.write("\n" + comb_line);
                y = yale.readLine();
                a = arg.readLine();

            }        
        }
        //write last lines
        while (y != null){
            out.write("\n" + y);
            y = yale.readLine();
        }
        while (a != null){
            out.write("\n" + a);
            a = arg.readLine();
        }
         out.close();
         yale.close();
         arg.close();
    }    
    public static void main(String[] args) throws FileNotFoundException, IOException{
        CombineYaleArgonne cya = new CombineYaleArgonne();
        //System.out.println(cya.checkSamplePositions("a\tb\tc\td\ty", "a\tb\tc\td"));
        //cya.checkSamplePositions("-	NA18486_yale	NA18498_yale	NA18499_yale	NA18501_yale	NA18502_yale	NA18504_yale	NA18505_yale	NA18507_yale	NA18508_yale	NA18510_yale	NA18511_yale	NA18516_yale	NA18517_yale	NA18519_yale	NA18520_yale	NA18522_yale	NA18523_yale	NA18852_yale	NA18853_yale	NA18855_yale	NA18856_yale	NA18858_yale	NA18859_yale	NA18861_yale	NA18862_yale	NA18870_yale	NA18871_yale	NA18909_yale	NA18912_yale	NA18913_yale	NA18916_yale	NA19093_yale	NA19098_yale	NA19099_yale	NA19101_yale	NA19102_yale	NA19108_yale	NA19114_yale	NA19116_yale	NA19119_yale	NA19127_yale	NA19128_yale	NA19130_yale	NA19131_yale	NA19137_yale	NA19138_yale	NA19140_yale	NA19141_yale	NA19143_yale	NA19144_yale	NA19147_yale	NA19152_yale	NA19153_yale	NA19159_yale	NA19160_yale	NA19171_yale	NA19172_yale	NA19190_yale	NA19192_yale	NA19193_yale	NA19200_yale	NA19201_yale	NA19203_yale	NA19204_yale	NA19207_yale	NA19209_yale	NA19210_yale	NA19222_yale	NA19225_yale	NA19238_yale", "-	NA18486_argonne	NA18498_argonne	NA18499_argonne	NA18501_argonne	NA18502_argonne	NA18504_argonne	NA18505_argonne	NA18507_argonne	NA18508_argonne	NA18510_argonne	NA18511_argonne	NA18516_argonne	NA18517_argonne	NA18519_argonne	NA18520_argonne	NA18522_argonne	NA18523_argonne	NA18852_argonne	NA18853_argonne	NA18855_argonne	NA18856_argonne	NA18858_argonne	NA18859_argonne	NA18861_argonne	NA18862_argonne	NA18870_argonne	NA18871_argonne	NA18909_argonne	NA18912_argonne	NA18913_argonne	NA18916_argonne	NA19093_argonne	NA19098_argonne	NA19099_argonne	NA19101_argonne	NA19102_argonne	NA19108_argonne	NA19114_argonne	NA19116_argonne	NA19119_argonne	NA19127_argonne	NA19128_argonne	NA19130_argonne	NA19131_argonne	NA19137_argonne	NA19138_argonne	NA19140_argonne	NA19141_argonne	NA19143_argonne	NA19144_argonne	NA19147_argonne	NA19152_argonne	NA19153_argonne	NA19159_argonne	NA19160_argonne	NA19171_argonne	NA19172_argonne	NA19190_argonne	NA19192_argonne	NA19193_argonne	NA19200_argonne	NA19201_argonne	NA19203_argonne	NA19204_argonne	NA19207_argonne	NA19209_argonne	NA19210_argonne	NA19222_argonne	NA19225_argonne	NA19238_argonne	NA19239_argonne	NA19257_argonne");
       /*
        cya.combine("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_coverageBed_exonwise_Pickrell_yale/expression_table.txt", 
                "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_coverageBed_exonwise_Pickrell_argonne/expression_table_cut.txt", 
                "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Yale+Argonne_exonwise/yale_argonne_expression.txt");
        */
        
        cya.combine("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Yale_noNorm/expression_table_Yale_noNorm_full.txt", 
                "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Argonne_noNorm/expression_table_Argonne_noNorm_full.txt", 
                "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Y+A_noNorm/expression_table_Y+A_noNorm_full.txt");
        
        
        //DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_coverageBed_Pickrell_yale/expression_table_full_cut.txt.10PCAsOverSamplesRemoved-GeneticVectorsNotRemoved_added.txt");
        //System.out.println(ds.colObjects.get(0));
        //System.out.println(ds.rowObjects.get(0));
                
    }
    
}
