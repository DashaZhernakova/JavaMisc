/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author dashazhernakova
 */
public class eQTLList {
    ArrayList<eQTL> eqtl_list;

    public eQTLList(ArrayList<eQTL> eqtlList) {
        eqtl_list = eqtlList;
    }
    public eQTLList(){
        
    }
    public ArrayList<eQTL> getWithoutLD(){
        //ArrayList<eQTL> new_list = eqtl_list;
        ArrayList<eQTL> new_list = new ArrayList<eQTL>();
        for (eQTL e : eqtl_list){
            if (! containsSNPinLD(new_list, e.gene, e.r_value) )
                new_list.add(e);
        }
        return new_list;
    }
    
    public ArrayList<eQTL> filterMostSignificantPerGene(){
        ArrayList<eQTL> new_list = new ArrayList<eQTL>();
        Map<String, eQTL> eqtls_map = new HashMap<String, eQTL>();
        for (eQTL e : eqtl_list){
            if (eqtls_map.containsKey(e.gene)){
                if (Math.abs(eqtls_map.get(e.gene).r_value) <  Math.abs(e.r_value))
                    eqtls_map.put(e.gene,e);
            }
            else
                eqtls_map.put(e.gene,e);
        }
        for (eQTL e : eqtls_map.values())
            new_list.add(e);
        return new_list;
    }
    public boolean containsSNPinLD(ArrayList<eQTL> eqtls, String gene, float r){
        for (eQTL e : eqtls){
            if ( (e.gene.equals(gene)) && (e.r_value == r) )
                return true;
        }
        return false;
    }
    
    public boolean contains(eQTL e, boolean probeLevel){
        boolean contains = false;
        
        if (probeLevel){
            for (eQTL eqtl : eqtl_list){
                if ( (eqtl.probe.equals(e.probe)) && (eqtl.snp.equals(e.snp)) )
                    return true;
            }
        }
        else{
            for (eQTL eqtl : eqtl_list){
                if ( (eqtl.gene.equals(e.gene)) && (eqtl.snp.equals(e.snp)) )
                    return true;
            }
        }
        
        return contains;
    }
    
    
}
