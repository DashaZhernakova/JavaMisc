/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package processtmap;
import java.io.File;
import java.io.IOException;
import umcg.genetica.io.text.TextFile;
/**
 *
 * @author dashazhernakova
 */
public class MakeMapPed {
    //ArrayList<String> ped;
    String[] ped;
    public MakeMapPed(){
        //ped = new ArrayList<String>();
        
    }
    
    public void generate(String in_f) throws IOException{
        TextFile in = new TextFile(in_f, false);
        TextFile ped_out = new TextFile(in_f.replace(".txt", ".ped"), true);
        TextFile map_out = new TextFile(in_f.replace(".txt", ".map"), true);
        in.open();
        ped_out.open();
        map_out.open();
        String[] els = in.readLineElems(TextFile.tab);
        int len = els.length;
        ped = new String[len];
        for (int i = 4; i < len; i++){
            //ped.add(els[i]); 
            ped[i] = Integer.toString(i - 3) + " " + els[i] + " 0 0 1 2";
        }
        while( (els = in.readLineElems(TextFile.tab)) != null ){            
            map_out.writeln(els[0] + " " + els[2] + " 0 " + els[1]);
            for (int i = 4; i < len; i ++){
                //String tmp = ped.get(i - 4);
                //tmp+=" "+els[i].charAt(0)+" "+els[i].charAt(1);
                ped[i] += " "+els[i].charAt(0)+" "+els[i].charAt(1);
            }
        }
        for (int i = 4; i< len; i ++){
            ped_out.writeln(ped[i]);
        }
        ped_out.close();
        map_out.close();
    }
    public void generateTransposedMontgomery(String in_f) throws IOException{
        TextFile in = new TextFile(in_f, false);
        TextFile ped_out = new TextFile(in_f.replace(".txt", ".ped"), true);
        TextFile map_out = new TextFile(in_f.replace(".txt", ".map"), true);
        in.open();
        ped_out.open();
        map_out.open();
        int startInfo = 4;
        String[] els = in.readLineElems(TextFile.tab);
        int len = els.length;
        ped = new String[len];
        ped_out.write("1");
        for (int i = 5; i < len; i++){
            //ped.add(els[i]); 
            //ped[i] = Integer.toString(i - 3) + " " + els[i] + " 0 0 1 2";
            ped_out.write(" " + Integer.toString(i - 3));
        }
        ped_out.write("\n" + els[startInfo]);
        for (int i = startInfo + 1; i < len; i++)
            ped_out.write(" " + els[i]);
        
        ped_out.write("\n0");
        for (int i = startInfo+1; i < len; i++)
            ped_out.write(" 0");
        
        ped_out.write("\n0");
        for (int i = startInfo+1; i < len; i++)
            ped_out.write(" 0"); 
        
        ped_out.write("\n1");
        for (int i = startInfo+1; i < len; i++)
            ped_out.write(" 1");
        
        ped_out.write("\n2");
        for (int i = startInfo+1; i < len; i++)
            ped_out.write(" 2");
        
        
        
        while( (els = in.readLineElems(TextFile.tab)) != null ){            
            map_out.writeln(els[0] + " " + els[2] + " 0 " + els[1]);
            String next = "";
            ped_out.write("\n" + String.valueOf(els[startInfo].charAt(0)));
            next+=String.valueOf(els[startInfo].charAt(1));
            for (int i = startInfo+1; i < len; i ++){
                //String tmp = ped.get(i - 4);
                //tmp+=" "+els[i].charAt(0)+" "+els[i].charAt(1);
                ped_out.write(" "+els[i].charAt(0));
                next+=" "+els[i].charAt(1);
            }
            ped_out.write("\n" + next);
        }
        
        ped_out.close();
        map_out.close();
    }
    public void generateTransposedPickrell(String dir_name) throws IOException{
        File dir = new File(dir_name);
        for (File ch : dir.listFiles()) {
            if (ch.getName().endsWith(".txt")){
                
                String in_f = ch.getPath();
                System.out.println("processing " + in_f);
                TextFile in = new TextFile(in_f, false);
                TextFile ped_out = new TextFile(in_f.replace(".txt", ".ped"), true);
                TextFile map_out = new TextFile(in_f.replace(".txt", ".map"), true);
                in.open();
                ped_out.open();
                map_out.open();
                int startInfo = 11;
                String[] els = in.readLineElems(TextFile.space);
                int len = els.length;
                ped = new String[len];
                ped_out.write("1");
                for (int i = startInfo+1; i < len; i++){
                    //ped.add(els[i]); 
                    //ped[i] = Integer.toString(i - 3) + " " + els[i] + " 0 0 1 2";
                    ped_out.write(" " + Integer.toString(i - 3));
                }
                ped_out.write("\n" + els[startInfo]);
                for (int i = startInfo+1; i < len; i++){
                    ped_out.write(" " + els[i]);
                }

                ped_out.write("\n0");
                for (int i = startInfo+1; i < len; i++)
                    ped_out.write(" 0");

                ped_out.write("\n0");
                for (int i = startInfo+1; i < len; i++)
                    ped_out.write(" 0"); 

                ped_out.write("\n1");
                for (int i = startInfo+1; i < len; i++)
                    ped_out.write(" 1");

                ped_out.write("\n2");
                for (int i = startInfo+1; i < len; i++)
                    ped_out.write(" 2");



                while( (els = in.readLineElems(TextFile.space)) != null ){            
                    map_out.writeln(els[2].replaceFirst("chr", "") + "\t" + els[0] + "\t0\t" + els[3]);
                    String next = "";
                    ped_out.write("\n" + String.valueOf(els[startInfo ].charAt(0)));
                    next+=String.valueOf(els[startInfo].charAt(1));
                    for (int i = startInfo+1; i < len; i ++){
                        //String tmp = ped.get(i - 4);
                        //tmp+=" "+els[i].charAt(0)+" "+els[i].charAt(1);
                        ped_out.write(" "+els[i].charAt(0));
                        next+=" "+els[i].charAt(1);
                    }
                    ped_out.write("\n" + next);
                }

                ped_out.close();
                map_out.close();
            }
        }
    }
    public void get(int[] posToLeave, String in_f, String out_f) throws IOException{
        TextFile in = new TextFile(in_f, false);
        TextFile out = new TextFile(out_f, true);
        in.open();
        out.open();
        String[] els;
        
        while( (els = in.readLineElems(TextFile.tab)) != null ){
            String new_str = "";
            for (int i = 0; i < posToLeave.length; i++)
                new_str += "\t" + els[posToLeave[i]];
            out.writeln(new_str.replaceFirst("\t", ""));
        }
        in.close();
        out.close();
    }
    

    
    public static void main(String[] args) throws IOException {
        MakeMapPed m = new MakeMapPed();
        int[] x = {2,3};
        m.generateTransposedPickrell("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/HapMap3YRI");
    }
}
