/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import com.lowagie.text.Anchor;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author dashazhernakova
 */
public class GeneNameConverter {
    Map<String, String> idsToHUGO = new HashMap<String, String>();
    Map<String, String> trIdsToGeneIds = new HashMap<String, String>();
    Map<String, String> HUGOToIds = new HashMap<String, String>();
    public GeneNameConverter() throws IOException {
        TextFile t = new TextFile("/Users/dashazhernakova/Documents/UMCG/hg19/ensembl_gene_names.txt", false);
        t.open();
        String[] els;
        
        while ( (els = t.readLineElems(TextFile.tab)) != null )
            idsToHUGO.put(els[0], els[1]);
        t.close();
    }
    
    public GeneNameConverter(String fName) throws IOException {
        TextFile t = new TextFile(fName, false);
        t.open();
        String[] els;
        
        while ( (els = t.readLineElems(TextFile.tab)) != null )
            idsToHUGO.put(els[0], els[1]);
        t.close();
    }
    public void makeIdToHUGOmap()throws IOException{
       TextFile t = new TextFile("/Users/dashazhernakova/Documents/UMCG/hg19/ensembl_gene_names.txt", false);
       String[] els;
       while ( (els = t.readLineElems(TextFile.tab)) != null )
            HUGOToIds.put(els[1], els[0]);
        t.close();
        
    }
    public void makeHugoToIdmap(){
        
    }
    public String getId(String hugo){
        return HUGOToIds.get(hugo);
    }
    public void toHUGO() throws IOException{
        makeIdToHUGOmap();
        TextFile tr = new TextFile("/Users/dashazhernakova/Documents/UMCG/hg19/ucsc_refGene_hg19.bed", false);
        TextFile out = new TextFile("/Users/dashazhernakova/Documents/UMCG/hg19/ucsc_refGene_HUGO_hg19.bed", true);
        String line = "";
        String[] spl;
        while ( (line = tr.readLine()) != null){
            spl = line.split("\t");
            spl[3] = idsToHUGO.get(spl[3]);
            out.writelnTabDelimited(spl);
        }
        out.close();
        tr.close();
    }
    
    public String getHugoName(String id){
        return idsToHUGO.get(id);
    }
    public ArrayList<String> geneIdsToTrIds(String gene) throws IOException{
        TextFile t = new TextFile("/Users/dashazhernakova/Documents/UMCG/hg19/transcrIds_to_geneIds.txt", false);
        t.open();
        String[] els;
        Map<String, ArrayList<String>> genesIdsToTrIds = new HashMap<String, ArrayList<String>>();
        while ( (els = t.readLineElems(TextFile.tab)) != null ){
            if (genesIdsToTrIds.containsKey(els[1]))
                genesIdsToTrIds.get(els[1]).add(els[0]);
            else{
                ArrayList l = new ArrayList<String>();
                l.add(els[0]);
                genesIdsToTrIds.put(els[1], l);
            }
                
        }
         t.close();
        return genesIdsToTrIds.get(gene);
       
    }
    public void transcrIdsToGeneNames() throws IOException{
        TextFile t = new TextFile("/Users/dashazhernakova/Documents/UMCG/hg19/transcrIds_to_geneIds.txt", false);
        t.open();
        String[] els;
        
        while ( (els = t.readLineElems(TextFile.tab)) != null )
            trIdsToGeneIds.put(els[0], els[1]);
        t.close();
        
        TextFile tr = new TextFile("/Users/dashazhernakova/Documents/UMCG/tmp_tr.txt", false);
        tr.open();
        TextFile out = new TextFile("/Users/dashazhernakova/Documents/UMCG/tmp_genes.txt", true);
        out.open();
        String line = "";
        while ( (line = tr.readLine()) != null){
            out.write(trIdsToGeneIds.get(line) + "\n");
        }
        out.close();
        tr.close();
    }
     public static void main(String[] args) throws IOException {
         GeneNameConverter c = new GeneNameConverter("/Users/dashazhernakova/Documents/UMCG/hg19/ucsc_refGeneIds2HUGO.txt");
         //c.transcrIdsToGeneNames();
         c.toHUGO();
         
     }
}
