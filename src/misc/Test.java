package misc;

import java.io.IOException;
import umcg.genetica.io.trityper.converters.TriTyperToVCF;

/**
 *
 * @author dashazhernakova
 */
public class Test {
    
    public void f(StringBuffer g){
        if (g != null)
            g.append(" world");
        else
            g = new StringBuffer( "Hello" );
        System.out.println("f: " + g);
    } 
    public void x(int i){
        i++;
        System.out.println("x: " + i);
    }
    
    public void swap(StringBuffer s1, StringBuffer s2){
        StringBuffer tmp;
        s1.append("xxx");
        tmp = s1;
        s1 = s2;
        s2 = tmp;
        System.out.println(s1 + ", " + s2);
        
    }
    public static void main(String[] args) throws IOException{
        TriTyperToVCF conv = new TriTyperToVCF();
        conv.convert("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/newGenotypes", "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/newGenotypes/", 
                "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/newGenotypes/chr6_100SNPs_list.txt");
    }
}
