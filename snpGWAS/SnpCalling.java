import java.io.*;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;

/**
 * Created by jessiewu on 15/07/15.
 */
public class SnpCalling {
    public static final String STATE_A = "A";
    public static final String STATE_C = "C";
    public static final String STATE_G = "G";
    public static final String STATE_T = "T";
    public static final String STATE_GAP = "-";
    public static final String STATE_N = "N";

    public static final String[] STATES = new String[]{STATE_A, STATE_C, STATE_G, STATE_T};

    public static void main(String[] args){
        String prefix = args[0];

        String reference = args[1];

        String dataFile = args[2];
        try{

            String refSeq = getSeqStringFromFasta(reference);
            int seqLength = refSeq.length();

            int[] countsA = new int[seqLength];
            int[] countsC = new int[seqLength];
            int[] countsG = new int[seqLength];
            int[] countsT = new int[seqLength];
            int[] countsGAP = new int[seqLength];
            int[] countsN = new int[seqLength];

            int[] presenceA = new int[seqLength];
            int[] presenceC = new int[seqLength];
            int[] presenceG = new int[seqLength];
            int[] presenceT = new int[seqLength];
            int[] presenceGAP = new int[seqLength];
            int[] presenceN = new int[seqLength];

            String[] filePaths = getFilePaths(dataFile);

            int seqCount = filePaths.length;



            for(int i = 0; i < seqCount;i++) {

                //System.out.println((0 > 5? 0:1));
                String filePath = filePaths[i];
                String seq = getSeqStringFromFastaGz(filePath);
                    //System.out.println(seq.length());

                    //String lineElt = everything.toString().split("");
                    //System.out.println(lineElt.length);

                    getCharIndexInSeq(seq, STATE_A, countsA);
                    getCharIndexInSeq(seq, STATE_C, countsC);
                    getCharIndexInSeq(seq, STATE_G, countsG);
                    getCharIndexInSeq(seq, STATE_T, countsT);
                    getCharIndexInSeq(seq, STATE_GAP, countsGAP);
                    getCharIndexInSeq(seq, STATE_N, countsN);




                    if(i % 100 == 0){
                        System.out.println(i);
                    }






            }

            int[] alleleTypeCount = getAlleleTypeCount(countsA, countsC, countsG, countsT);

            Integer[][] alleleTypeIndex = getBiallelicSnpIndices(alleleTypeCount);
            Integer[] polyIndex = alleleTypeIndex[0];
            Integer[] biallelicIndex =  alleleTypeIndex[1];
            Integer[] tritetrallelicIndex =  alleleTypeIndex[2];

            System.out.println("Number of biallelic sites: "+biallelicIndex.length);
            System.out.println("Number of tritetrallelic sites: "+ tritetrallelicIndex.length);





            //String[][] alleleId  = new String[polyIndex.length][4];

            /*for(int iPoly = 0; iPoly < polyIndex.length; iPoly++){

                 alleleId[iPoly] = sortState(countsA[polyIndex[iPoly]],
                         countsC[polyIndex[iPoly]],
                         countsG[polyIndex[iPoly]],
                         countsT[polyIndex[iPoly]]
                 );



             }*/

        String[][] biallelicState  = new String[biallelicIndex.length][4];
        String[][] tritetrallelicState  = new String[tritetrallelicIndex.length][4];

        for(int i = 0; i < biallelicIndex.length; i++){

            biallelicState[i] = sortState(countsA[biallelicIndex[i]]+0.3,
                    countsC[biallelicIndex[i]]+0.2,
                    countsG[biallelicIndex[i]]+0.1,
                    countsT[biallelicIndex[i]]
            );
        }

        for(int i = 0; i < tritetrallelicIndex.length; i++){

            tritetrallelicState[i] = sortState(countsA[tritetrallelicIndex[i]]+0.3,
                    countsC[tritetrallelicIndex[i]]+0.2,
                    countsG[tritetrallelicIndex[i]]+0.1,
                    countsT[tritetrallelicIndex[i]]
            );
        }

        int[][] bip = new int[seqCount][biallelicIndex.length];
        int[][] snp = new int[seqCount][tritetrallelicIndex.length];


            for(int iSeq = 0; iSeq < seqCount; iSeq++){

                String filePath = filePaths[iSeq];
                String seq = getSeqStringFromFastaGz(filePath);


                bip[iSeq] = getBiallelicSnpPatterns(seq, biallelicIndex, biallelicState);
                snp[iSeq] = getTritetrallelicSnpPatterns(seq, tritetrallelicIndex, tritetrallelicState);

                if(iSeq % 100 == 0){
                    System.out.println("Stage2: "+iSeq);
                }


            }

            writeInfo(prefix+".bipinfo.txt",
                    prefix+".ttpinfo.txt",
                    biallelicIndex,
                    biallelicState,
                    tritetrallelicIndex,
                    tritetrallelicState,
                    countsA, countsC, countsG, countsT
                    );

            String[] ids = getIds(dataFile);

            writePatterns(prefix+".bip.patterns.txt",
                    prefix+".ttp.patterns.txt",
                    bip,
                    snp,
                    ids
            );

        }catch(Exception e){
            throw new RuntimeException(e);
        }

        System.out.println("Completed");


    }

    public static String[] getFilePaths(String dataFilePath){
        ArrayList<String> filePathList = new ArrayList<String>();

        try{
            BufferedReader reader = new BufferedReader(new FileReader(dataFilePath));
            String line = reader.readLine();
            int index = Arrays.asList(line.split("\t")).indexOf("filePath");

            while((line = reader.readLine()) != null){
                filePathList.add(line.trim().split("\t")[index]);
            }

            reader.close();

        }catch(Exception e ){
            throw new RuntimeException(e);
        }

        String[] filePaths = new String[filePathList.size()];

        filePaths = filePathList.toArray(filePaths);

        return filePaths;

    }


    public static String[] getIds(String dataFilePath){
        ArrayList<String> idList = new ArrayList<String>();

        try{
            BufferedReader reader = new BufferedReader(new FileReader(dataFilePath));
            String line = reader.readLine();
            int index = Arrays.asList(line.split("\t")).indexOf("id");

            while((line = reader.readLine()) != null){
                idList.add(line.trim().split("\t")[index]);
            }

            reader.close();

        }catch(Exception e ){
            throw new RuntimeException(e);
        }

        String[] ids = new String[idList.size()];

        ids = idList.toArray(ids);

        return ids;

    }


    public static String getSeqStringFromFastaGz(String filePath){
        StringBuilder everything = new StringBuilder();

        try {

            FileInputStream fin = new FileInputStream(filePath);
            GZIPInputStream gzis = new GZIPInputStream(fin);
            InputStreamReader xover = new InputStreamReader(gzis);
            BufferedReader seqReader = new BufferedReader(xover);

            String line = seqReader.readLine();
            while ((line = seqReader.readLine()) != null) {
                everything.append(line);
            }

            fin.close();
            gzis.close();
            xover.close();
            seqReader.close();


        }catch(Exception e){
            throw new RuntimeException(e);
        }

        return everything.toString();

    }

    public static String getSeqStringFromFasta(String filePath){
        System.out.println(filePath);
        StringBuilder everything = new StringBuilder();

        try {

            BufferedReader seqReader = new BufferedReader(new FileReader(filePath));

            String line = seqReader.readLine();
            while ((line = seqReader.readLine()) != null) {
                everything.append(line);
            }


            seqReader.close();


        }catch(Exception e){
            throw new RuntimeException(e);
        }

        return everything.toString();

    }


    public static void getCharIndexInSeq(String seq, String state, int[] counts){
        for (int index = seq.indexOf(state); index >= 0; index = seq.indexOf(state, index + 1)) {
            if(index > -1) {
                counts[index]++;
                //System.out.println(index);
            }
        }

    }

    public static String[] sortState(double countA, double countC, double countG, double countT){
        double[] counts = new double[]{countA, countC, countG, countT};
        HashMap<Double,Integer> countMap = new HashMap<Double,Integer>();
        countMap.put(countA, 0);
        countMap.put(countC, 1);
        countMap.put(countG, 2);
        countMap.put(countT, 3);
        Arrays.sort(counts);
        String[] sortedStates = new String[4];
        sortedStates[0] = STATES[countMap.get(counts[3])];
        sortedStates[1] = STATES[countMap.get(counts[2])];
        sortedStates[2] = STATES[countMap.get(counts[1])];
        sortedStates[3] = STATES[countMap.get(counts[0])];

        return sortedStates;

    }

    public static int[]  getAlleleTypeCount(int[] countsA, int[] countsC, int[] countsG, int[] countsT){
        int[] alleleTypeCounts = new int[countsA.length];
        int presA, presC, presG, presT;
        for(int i = 0; i < alleleTypeCounts.length; i++){

            presA = countsA[i] > 0? 1:0;
            presC = countsC[i] > 0? 1:0;
            presG = countsG[i] > 0? 1:0;
            presT = countsT[i] > 0? 1:0;
            alleleTypeCounts[i] = presA + presC + presG + presT;



        }

        return alleleTypeCounts;



    }

    public static Integer[][] getBiallelicSnpIndices(int[] alleleTypeCount){
        ArrayList<Integer> isPoly = new ArrayList<Integer>();
        ArrayList<Integer> isBiallelic = new ArrayList<Integer>();
        ArrayList<Integer> isTriTetrallic = new ArrayList<Integer>();

        for(int i = 0; i < alleleTypeCount.length; i++){
            if(alleleTypeCount[i] == 2){
                isPoly.add(i);
                isBiallelic.add(i);
            }else if(alleleTypeCount[i] == 3 || alleleTypeCount[i] == 4){
                isPoly.add(i);
                isTriTetrallic.add(i);
            }else if(alleleTypeCount[i] >= 5 || alleleTypeCount[i] < 0){
                throw new RuntimeException("Allele type count is "+alleleTypeCount[i]+" at index "+i+". ");
            }
        }
        Integer[] polyIndex = isPoly.toArray(new Integer[isPoly.size()]);
        Integer[] biallelicIndex = isBiallelic.toArray(new Integer[isBiallelic.size()]);
        Integer[] tritetraIndex = isTriTetrallic.toArray(new Integer[isTriTetrallic.size()]);
        Integer[][] alleleTypeIndex = new Integer[3][];
        alleleTypeIndex[0] = polyIndex;
        alleleTypeIndex[1] = biallelicIndex;
        alleleTypeIndex[2] = tritetraIndex;



        return alleleTypeIndex;


    }
    /*# Second pass over the files: populate bip and snp objects
    for(i in 1:length(filepaths)) {
        message(filepaths[i])
        # Read the mapcall file
        fa = toupper(read_reference_gz(filepaths[i]))
        fa.poly = fa[is.poly]
        # Populate bip object
        bip[,i] = -1
        bip[fa[isbiallelic]==allele.id[1,isbiallelic[is.poly]],i] = 0
        bip[fa[isbiallelic]==allele.id[2,isbiallelic[is.poly]],i] = 1
        # Populate snp object
        snp[,i] = -1
        snp[fa[is.poly & !isbiallelic]==allele.id[1,!isbiallelic[is.poly]],i] = 0
        snp[fa[is.poly & !isbiallelic]==allele.id[2,!isbiallelic[is.poly]],i] = 1
        snp[fa[is.poly & !isbiallelic]==allele.id[3,!isbiallelic[is.poly]],i] = 2
        snp[fa[is.poly & !isbiallelic]==allele.id[4,!isbiallelic[is.poly]],i] = 3
        cat("Done",i,"\n")
    }*/

    public static int[] getBiallelicSnpPatterns(String seq,
                                      Integer[] biallelicIndex,
                                      String[][] biallelicState){

        int[] pattern = new int[biallelicIndex.length];

        for(int i = 0; i < biallelicIndex.length; i++){
            if((""+seq.charAt(biallelicIndex[i])).equals(biallelicState[i][0])){
                pattern[i] = 0;

            }else if((""+seq.charAt(biallelicIndex[i])).equals(biallelicState[i][1])){
                pattern[i] = 1;

            }else{
                pattern[i] = -1;

            }
        }
        return pattern;

    }


    public static int[] getTritetrallelicSnpPatterns(String seq,
                                               Integer[] triterallelicIndex,
                                               String[][] triterallelicState){

        int[] pattern = new int[triterallelicIndex.length];

        for(int i = 0; i < triterallelicIndex.length; i++){
            if((""+seq.charAt(triterallelicIndex[i])).equals(triterallelicState[i][0])){

                pattern[i] = 0;

            }else if((""+seq.charAt(triterallelicIndex[i])).equals(triterallelicState[i][1])) {

                pattern[i] = 1;

            }else if((""+seq.charAt(triterallelicIndex[i])).equals(triterallelicState[i][2])) {

                pattern[i] = 2;

            }else if((""+seq.charAt(triterallelicIndex[i])).equals(triterallelicState[i][3])){

                pattern[i] = 3;

            }else{
                pattern[i] = -1;

            }
        }

        return pattern;

    }

    public static void writeInfo(String bipinfoOutputFilePath,
                                    String snpinfoOutputFilePath,
                                    Integer[] biallelicIndex,
                                    String[][] biallelicState,
                                    Integer[] tritetrallelicIndex,
                                    String[][] tritetrallelicState,
                                    int[] countsA,
                                    int[] countsC,
                                    int[] countsG,
                                    int[] countsT){
        try {

            PrintWriter writer = new PrintWriter(bipinfoOutputFilePath);
            writer.println("Position\t"+"Allele0\t"+"Allele1\t"+"A\t"+"C\t"+"G\t"+"T");
            for(int i = 0; i < biallelicIndex.length; i++){
                writer.println((biallelicIndex[i]+1) +"\t" +
                        biallelicState[i][0]+ "\t" +
                        biallelicState[i][1]+"\t" +
                        countsA[biallelicIndex[i]] + "\t" +
                        countsC[biallelicIndex[i]] + "\t" +
                        countsG[biallelicIndex[i]] + "\t" +
                        countsT[biallelicIndex[i]]
                );
            }

            writer.close();


            writer = new PrintWriter(snpinfoOutputFilePath);
            writer.println("Position\t"+"Allele0\t"+"Allele1\t"+"Allele2\t"+"Allele3\t"+"A\t"+"C\t"+"G\t"+"T");
            for(int i = 0; i < tritetrallelicIndex.length; i++){
                writer.println((tritetrallelicIndex[i]+1) +"\t" +
                                tritetrallelicState[i][0]+ "\t" +
                                tritetrallelicState[i][1]+"\t" +
                                tritetrallelicState[i][2]+"\t" +
                                tritetrallelicState[i][3]+"\t" +
                                countsA[tritetrallelicIndex[i]] + "\t" +
                                countsC[tritetrallelicIndex[i]] + "\t" +
                                countsG[tritetrallelicIndex[i]] + "\t" +
                                countsT[tritetrallelicIndex[i]]
                );
            }

            writer.close();


        }catch(Exception e){
            throw new RuntimeException(e);
        }





    }

    public static void writePatterns(String bippatternsPath,
                                     String snppatternsPath,
                                     int[][] bip,
                                     int[][]snp,
                                     String id[]){
        try{
            PrintWriter writer = new PrintWriter(bippatternsPath);
            for(int i = 0; i < id.length; i++){
                writer.print(id[i]+"\t");

            }
            writer.println();

            int bipCount = bip[0].length;

            for(int i = 0; i < bipCount; i++ ){
                for(int j = 0; j < bip.length; j++){
                    writer.print(bip[j][i]+"\t");
                }
                writer.println();
            }

            writer.close();


            writer = new PrintWriter(snppatternsPath);
            for(int i = 0; i < id.length; i++){
                writer.print(id[i]+"\t");

            }
            writer.println();

            int snpCount = snp[0].length;

            for(int i = 0; i < snpCount; i++ ){

                for(int j = 0; j < snp.length; j++){
                    writer.print(snp[j][i]+"\t");
                }
                writer.println();
            }
            writer.close();

        }catch(Exception e){
            throw new RuntimeException(e);
        }
    }







}
