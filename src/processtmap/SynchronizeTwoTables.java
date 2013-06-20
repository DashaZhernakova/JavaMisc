package processtmap;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author dashazhernakova
 */
public class SynchronizeTwoTables {
    //int nCol;
    public SynchronizeTwoTables(){
        //nCol = numCol;
    }
    public void synchronize(String f1, String f2, String outF) throws IOException{
        TextFile file1 = new TextFile(f1, false);
        TextFile file2 = new TextFile(f2, false);
        TextFile out = new TextFile(outF, true);
        Map<String,String> first = new HashMap<String, String>();
        
        String line;
        while ( (line = file1.readLine()) != null ){
            first.put(line.split("\t")[0], line);
        }
        
        while ((line = file2.readLine()) != null ){
            if (first.containsKey(line.split("\t")[0]))
                out.writeln(line + "\t" + first.get(line.split("\t")[0]));
        }
        
        file1.close();
        file2.close();
        out.close();
    }
    public static void main(String[] args) throws IOException {
        SynchronizeTwoTables s = new SynchronizeTwoTables();
        s.synchronize("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp1.txt", 
                "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp2.txt", 
                "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp_joint.txt");
    }
}
