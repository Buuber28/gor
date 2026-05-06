import java.io.*;
import java.util.*;

public class GOR4 {

    static final int window = 17;
    static final int half = window / 2;


    //reads a multiple alignment file into a list of FastaEntry objects
    //frist entry is the target and the following ones are the homologs
    private static List<FastaEntry> readMaf(String maPath) {
        List<FastaEntry> out = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(maPath))) {

            String line;

            String targetId = null;
            String targetSeq = null;

            while ((line = br.readLine()) != null) {

                line = line.trim();
                if (line.isEmpty()) continue;

                if (line.startsWith(">")) {
                    targetId = line.substring(1).trim();
                    continue;
                }


                if (line.startsWith("AS")) {

                    targetSeq = line.substring(2).trim().toUpperCase();
                    continue;
                }

                // just for testing with the gor_example files
                if (line.startsWith("SS")) {
                    continue;
                }

                String[] tok = line.split("\\s+");
                if (tok.length == 0) continue;

                String homId;
                String homSeq;

                if (tok.length >= 2 && tok[0].matches("\\d+")) {
                    homId = tok[0];
                    homSeq = tok[1];
                } else {
                    homId = String.valueOf(out.size() + 1);
                    homSeq = tok[0];
                }

                homSeq = homSeq.trim().toUpperCase();

                // if there is no sequence but a number
                if (homSeq.isEmpty() || homSeq.matches("\\d+")) continue;
                out.add(new FastaEntry(homId, homSeq));

            }

            if (targetId == null || targetSeq == null) {
                throw new RuntimeException("Bad alignment file");
            }

            out.add(0, new FastaEntry(targetId, targetSeq));


        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        return out;

    }

    //read fasta sequence and save it as FastaEntry
    private static List<FastaEntry> readFasta(String fastaPath) {
        List<FastaEntry> out = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(fastaPath))) {
            String line;
            String currentId = null;
            StringBuilder currentSeq = new StringBuilder();

            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                if (line.startsWith(">")) {
                    // flush previous
                    if (currentId != null) {
                        out.add(new FastaEntry(currentId, currentSeq.toString()));
                    }
                    currentId = line.substring(1).trim();
                    currentSeq.setLength(0);
                } else {
                    currentSeq.append(line);
                }
            }

            // flush last
            if (currentId != null) {
                out.add(new FastaEntry(currentId, currentSeq.toString()));
            }

        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        return out;
    }

    //reads model file and fills Gor4modelData with raw counts for 6D and 4D matrices
    public static Gor4ModelData readModel(String modelPath) {
        final String AA = "ACDEFGHIKLMNPQRSTVWY";
        Map<Character, Integer> aaMap = aa_number(AA);

        Gor4ModelData data = new Gor4ModelData(AA);

        try (BufferedReader br = new BufferedReader(new FileReader(modelPath))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;
                if (line.startsWith("//")) continue;

                if (line.startsWith("=") && line.endsWith("=")) {
                    String header = line.substring(1, line.length() - 1).trim();
                    String[] parts = header.split(",");

                    // 6D: =State,CenterAA,AAP,pOff=
                    if (parts.length == 4) {
                        char stateChar = parts[0].trim().charAt(0);
                        char centerAAChar = parts[1].trim().charAt(0);
                        char aaPChar = parts[2].trim().charAt(0);
                        int pOff = Integer.parseInt(parts[3].trim());
                        int pIdx = pOff + 8;

                        int s = stateIndex(stateChar);
                        Integer centerAA = aaMap.get(centerAAChar);
                        Integer aaP = aaMap.get(aaPChar);

                        if (s < 0 || centerAA == null || aaP == null || pIdx < 0 || pIdx > 16) {
                            skipTable(br, 20);
                            continue;
                        }

                        int rowsRead = 0;
                        while (rowsRead < 20) {
                            String row = br.readLine();
                            if (row == null) throw new RuntimeException("Unexpected EOF in 6D table");
                            row = row.trim();
                            if (row.isEmpty() || row.startsWith("//")) continue;

                            String[] tok = row.split("\\s+");
                            if (tok.length < 18) throw new RuntimeException("Bad 6D row: " + row);

                            char aaQChar = tok[0].charAt(0);
                            Integer aaQ = aaMap.get(aaQChar);
                            if (aaQ == null) throw new RuntimeException("Unknown AA in 6D row: " + aaQChar);

                            for (int qIdx = 0; qIdx < 17; qIdx++) {
                                long v = Long.parseLong(tok[1 + qIdx]);
                                data.counts6D[s][centerAA][aaP][pIdx][aaQ][qIdx] = v;
                            }
                            rowsRead++;
                        }
                    }

                    // 4D: =CenterAA,State=
                    else if (parts.length == 2) {
                        char centerAAChar = parts[0].trim().charAt(0);
                        char stateChar = parts[1].trim().charAt(0);

                        int s = stateIndex(stateChar);
                        Integer centerAA = aaMap.get(centerAAChar);

                        if (s < 0 || centerAA == null) {
                            skipTable(br, 20);
                            continue;
                        }

                        int rowsRead = 0;
                        while (rowsRead < 20) {
                            String row = br.readLine();
                            if (row == null) throw new RuntimeException("Unexpected EOF in 4D table");
                            row = row.trim();
                            if (row.isEmpty() || row.startsWith("//")) continue;

                            String[] tok = row.split("\\s+");
                            if (tok.length < 18) throw new RuntimeException("Bad 4D row: " + row);

                            char aaPChar = tok[0].charAt(0);
                            Integer aaP = aaMap.get(aaPChar);
                            if (aaP == null) throw new RuntimeException("Unknown AA in 4D row: " + aaPChar);

                            for (int pIdx = 0; pIdx < 17; pIdx++) {
                                long v = Long.parseLong(tok[1 + pIdx]);
                                data.counts4D[s][centerAA][aaP][pIdx] = v;
                            }
                            rowsRead++;
                        }
                    } else {
                        throw new RuntimeException("Unknown header format: " + line);
                    }
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        return data;
    }

    // reads the .db and turns entries into SeclibEntry Objects
    public static List<SeclibEntry> readSequences(String seclibFile) {
        List<SeclibEntry> entries = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(seclibFile))) {
            String line;
            String currentId = null;
            String currentAS = null;
            String currentSS = null;

            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                if (line.startsWith(">")) {
                    // flush previous entry
                    if (currentId != null) {

                        entries.add(new SeclibEntry(currentId, currentAS, currentSS));
                    }
                    currentId = line.substring(1).trim();
                    currentAS = null;
                    currentSS = null;

                } else if (line.startsWith("AS")) {
                    currentAS = line.substring(2).trim();

                } else if (line.startsWith("SS")) {
                    currentSS = line.substring(2).trim();

                }
            }

            // flush last entry
            if (currentId != null) {
                entries.add(new SeclibEntry(currentId, currentAS, currentSS));
            }

        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        return entries;
    }

    // write the matrix (raw counts)
    // counts6D[state][centerAA][aaP][pIdx][aaQ][qIdx]
    public static void writeModel6D(String modelFile, String aminoAcids, long[][][][][][] counts6D) {
        try (PrintWriter w = new PrintWriter(new FileWriter(modelFile))) {
            w.println("// Matrix6D");
            w.println();

            char[] states = new char[]{'C', 'E', 'H'};

            for (int s = 0; s < 3; s++) {
                for (int centerAA = 0; centerAA < 20; centerAA++) {
                    for (int aaP = 0; aaP < 20; aaP++) {
                        for (int pIdx = 0; pIdx < window; pIdx++) {

                            int pOff = pIdx - half;
                            w.println("=" + states[s] + "," +
                                    aminoAcids.charAt(centerAA) + "," +
                                    aminoAcids.charAt(aaP) + "," +
                                    pOff + "=");
                            w.println();

                            for (int aaQ = 0; aaQ < 20; aaQ++) {
                                w.print(aminoAcids.charAt(aaQ));
                                w.print("\t");
                                for (int qIdx = 0; qIdx < window; qIdx++) {
                                    w.print(counts6D[s][centerAA][aaP][pIdx][aaQ][qIdx]);
                                    w.print("\t");
                                }
                                w.println();
                            }
                            w.println();
                        }
                    }
                }
            }

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


    // write the 4D matrix (raw counts)
    // counts4D[state][centerAA][aaP][pIdx]
    // header: =<centerAA>,<state>=
    public static void writeModel4D(String modelFile, String aminoAcids, long[][][][] counts4D) {
        try (PrintWriter w = new PrintWriter(new FileWriter(modelFile, true))) {

            w.println("// Matrix4D");
            w.println();

            char[] states = new char[]{'C', 'E', 'H'};

            for (int centerAA = 0; centerAA < 20; centerAA++) {
                for (int s = 0; s < 3; s++) {

                    w.println("=" + aminoAcids.charAt(centerAA) + "," + states[s] + "=");
                    w.println();

                    for (int aaP = 0; aaP < 20; aaP++) {
                        w.print(aminoAcids.charAt(aaP));
                        w.print("\t");
                        for (int pIdx = 0; pIdx < window; pIdx++) {
                            w.print(counts4D[s][centerAA][aaP][pIdx]);
                            w.print("\t");
                        }
                        w.println();
                    }

                    w.println();
                }
            }

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


    private static void runTrain(TrainArgs t) {
        String aminoAcids = "ACDEFGHIKLMNPQRSTVWY";
        Map<Character, Integer> aaIndex = aa_number(aminoAcids);

        List<SeclibEntry> entries = readSequences(t.dbPath);

        Gor4Trainer trainer = new Gor4Trainer(aaIndex);
        Gor4Model model = trainer.train(entries);

        writeModel6D(t.modelPath, aminoAcids, model.getCounts6D());
        writeModel4D(t.modelPath, aminoAcids, model.getCounts4D());

    }

    private static void runPredict(PredictArgs p) {

        // read raw counts from model file
        Gor4ModelData data = readModel(p.modelPath);

        // build model + load counts + compute weights
        Map<Character, Integer> aaIndex = aa_number(data.aaOrder);

        Gor4Model model = new Gor4Model(half, aaIndex);

        model.loadCounts6D(data.counts6D);
        model.loadCounts4D(data.counts4D);
        model.finalizeModel();

        // read fasta and predict
        List<FastaEntry> fasta = readFasta(p.fastaPath);

        for (FastaEntry fe : fasta) {
            String id = fe.id;
            String seq = fe.seq.trim().toUpperCase();

            String ss = model.predictPS(seq);
            ss = postProcess(ss);

            if (ss.length() != seq.length()) {
                System.err.println("LENGTH MISMATCH id=" + id);
                System.err.println("AS len=" + seq.length());
                System.err.println("PS len=" + ss.length());
                System.err.println("AS=" + seq);
                System.err.println("PS=" + ss);
                throw new RuntimeException("PS length != AS length for " + id);
            }


            System.out.println(">" + id);
            System.out.println("AS " + seq);
            System.out.println("PS " + ss);


            if (p.probabilities) {

                // Output digits 0-9 scaled from probabilities
                StringBuilder pH = new StringBuilder(seq.length());
                StringBuilder pE = new StringBuilder(seq.length());
                StringBuilder pC = new StringBuilder(seq.length());

                for (int j = 0; j < seq.length(); j++) {
                    if (j < half || j >= seq.length() - half) {
                        pH.append('-');
                        pE.append('-');
                        pC.append('-');
                        continue;
                    }

                    double[] sc = model.scoreAt(seq, j); //

                    // softmax
                    double m = Math.max(sc[0], Math.max(sc[1], sc[2]));
                    double e0 = Math.exp(sc[0] - m);
                    double e1 = Math.exp(sc[1] - m);
                    double e2 = Math.exp(sc[2] - m);
                    double sum = e0 + e1 + e2;

                    double pCval = e0 / sum;
                    double pEval = e1 / sum;
                    double pHval = e2 / sum;

                    pC.append(probToDigit(pCval));
                    pE.append(probToDigit(pEval));
                    pH.append(probToDigit(pHval));
                }

                System.out.println("PH " + pH);
                System.out.println("PE " + pE);
                System.out.println("PC " + pC);
            }
        }
    }



    private static void runPredictMaf(PredictArgs p) {

        // load model counts and build weights
        Gor4ModelData data = readModel(p.modelPath);
        Map<Character, Integer> aaIndex = aa_number(data.aaOrder);

        Gor4Model model = new Gor4Model(half, aaIndex);
        model.loadCounts4D(data.counts4D);
        model.loadCounts6D(data.counts6D);
        model.finalizeModel();

        // list alignment files
        File dir = new File(p.mafPath);
        File[] files = dir.listFiles();
        if (files == null) throw new RuntimeException("Could not list files in --maf: " + p.mafPath);
        Arrays.sort(files);

        for (File f : files) {
            if (!f.isFile()) continue;

            String name = f.getName().toLowerCase();
            if (!(name.endsWith(".aln") || name.endsWith(".ma") || name.endsWith(".maf"))) continue;

            List<FastaEntry> aln = readMaf(f.getAbsolutePath());
            FastaEntry target = aln.get(0);
            String targetSeq = target.seq;
            int L = targetSeq.length();

            for (FastaEntry e : aln) {
                if (e.seq.length() != L) {
                    throw new RuntimeException("Alignment length mismatch in " + f.getName()
                            + ": target=" + L + " but " + e.id + "=" + e.seq.length());
                }
            }



            StringBuilder AS = new StringBuilder();
            StringBuilder PS = new StringBuilder();

            StringBuilder PH = new StringBuilder();
            StringBuilder PE = new StringBuilder();
            StringBuilder PC = new StringBuilder();

            // i
            for (int j = 0; j < L; j++) {
                char t = targetSeq.charAt(j);

                // skip gaps
                if (t == '-') continue;

                AS.append(t);

                if (!aaIndex.containsKey(t)) {
                    PS.append('C');
                    if (p.probabilities) { PH.append('0'); PE.append('0'); PC.append('9'); }
                    continue;
                }


                if (j < half || j >= L - half) {
                    PS.append('-');
                    if (p.probabilities) { PH.append('-'); PE.append('-'); PC.append('-'); }
                    continue;
                }

                // average P across sequences
                double sumPC = 0.0, sumPE = 0.0, sumPH = 0.0;
                int used = 0;

                for (FastaEntry e : aln) {
                    char c = e.seq.charAt(j);
                    if (c == '-') continue;                 // ignore gaps
                    if (!aaIndex.containsKey(c)) continue;  // ignore unknowns in homologs too

                    double[] sc = model.scoreAt(e.seq, j);

                    // softmax per seq -> probabilities
                    double m = Math.max(sc[0], Math.max(sc[1], sc[2]));
                    double e0 = Math.exp(sc[0] - m);
                    double e1 = Math.exp(sc[1] - m);
                    double e2 = Math.exp(sc[2] - m);
                    double all = e0 + e1 + e2;

                    double pCval = e0 / all;
                    double pEval = e1 / all;
                    double pHval = e2 / all;

                    sumPC += pCval;
                    sumPE += pEval;
                    sumPH += pHval;
                    used++;
                }

                if (used == 0) {
                    PS.append('C');
                    if (p.probabilities) { PH.append('0'); PE.append('0'); PC.append('9'); }
                    continue;
                }

                // average probs in this column
                double avgPC = sumPC / used;
                double avgPE = sumPE / used;
                double avgPH = sumPH / used;

                // pick max probability
                if (avgPH >= avgPE && avgPH >= avgPC) PS.append('H');
                else if (avgPE >= avgPH && avgPE >= avgPC) PS.append('E');
                else PS.append('C');

                if (p.probabilities) {
                    PC.append(probToDigit(avgPC));
                    PE.append(probToDigit(avgPE));
                    PH.append(probToDigit(avgPH));
                }



            }

            // overriding for window edges
            for (int k = 0; k < half && k < PS.length(); k++) {
                PS.setCharAt(k, '-');
                PS.setCharAt(PS.length() - 1 - k, '-');

                if (p.probabilities) {
                    PH.setCharAt(k, '-');
                    PE.setCharAt(k, '-');
                    PC.setCharAt(k, '-');

                    PH.setCharAt(PH.length() - 1 - k, '-');
                    PE.setCharAt(PE.length() - 1 - k, '-');
                    PC.setCharAt(PC.length() - 1 - k, '-');
                }
            }

            String psFinal = postProcess(PS.toString());


            // output
            System.out.println(">" + target.id);
            System.out.println("AS " + AS);
            System.out.println("PS " + psFinal);

            if (p.probabilities) {
                System.out.println("PH " + PH);
                System.out.println("PE " + PE);
                System.out.println("PC " + PC);
            }
        }
    }


    public static void main(String[] args) {

        // if training mode
        if (hasArg(args, "--db")) {
            TrainArgs t = parseTrainArgs(args);
            runTrain(t);
            return;
        }

        // otherwise prediction mode
        PredictArgs p = parsePredictArgs(args);


        if (p.modelPath == null || p.format == null) {
            printHelp();
            return;
        }

        // Exactly one of --seq or --maf mac
        boolean hasSeq = (p.fastaPath != null);
        boolean hasMaf = (p.mafPath != null);

        if (hasSeq == hasMaf) { // both true or both false
            System.err.println("Error: only exactly one of --seq or --maf allowed");
            printHelp();
            return;
        }

        if (hasSeq) {
            runPredict(p);
        } else {
            runPredictMaf(p);       // GOR5
        }
    }

    private static String postProcess(String ps){
        char[]  a = ps.toCharArray();


        // use a minimun length for H >= 3, E>=2
        a = enforceMin(a, 'H', 3);
        a = enforceMin(a, 'E',2);


        //fill single-gap inside long same state segments

        for(int i = 1; i < a.length -1; i++){

            if (a[i] == 'C' && a[i - 1] == a[i + 1] && (a[i - 1] == 'H' || a[i - 1] == 'E')) {
                a[i] = a[i-1];
            }
        }


        return new String(a);


    }

    private static char[] enforceMin(char[] a, char state, int min){
        int n = a.length;
        int i = 0;
        while(i<n){
            if(a[i] != state){
                i++ ;
                continue;
            }
            int j = i;
            while(j < n && a[j] == state) {j++;}
            int len = j - i;
            if(len<min){
                for (int k = i; k <j; k++) {a[k] = 'C';}
            }
            i = j;
        }
        return a;
    }


    //HELPERS

    public static int stateIndex(char s) {
        return switch (s) {
            case 'C' -> 0;
            case 'E' -> 1;
            case 'H' -> 2;
            default -> -1;
        };
    }

    //skips over the table, used when we find bad headers
    public static void skipTable(BufferedReader br, int rows) throws IOException {
        int read = 0;
        while (read < rows) {
            String row = br.readLine();
            if (row == null) return;
            row = row.trim();
            if (row.isEmpty() || row.startsWith("//")) continue;
            read++;
        }
    }

    // maps a AA to a number based on aa string order
    public static Map<Character, Integer> aa_number(String aminoAcids) {
        Map<Character, Integer> aaIndex = new HashMap<>();
        for (int i = 0; i < aminoAcids.length(); i++) {
            aaIndex.put(aminoAcids.charAt(i), i);
        }
        return aaIndex;
    }

    private static TrainArgs parseTrainArgs(String[] args) {
        TrainArgs t = new TrainArgs();

        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "--db" -> {
                    if (i + 1 >= args.length) {
                        printHelp();
                        return null;
                    }
                    t.dbPath = args[++i];
                }
                case "--method" -> {
                    if (i + 1 >= args.length) {
                        printHelp();
                        return null;
                    }
                    t.method = args[++i];
                }
                case "--model" -> {
                    if (i + 1 >= args.length) {
                        printHelp();
                        return null;
                    }
                    t.modelPath = args[++i];
                }
                default -> {

                }
            }
        }

        if (t.dbPath == null || t.method == null || t.modelPath == null) {
            printHelp();
            return null;
        }

        return t;
    }

    private static PredictArgs parsePredictArgs(String[] args) {
        PredictArgs p = new PredictArgs();

        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "--model" -> {
                    if (i + 1 >= args.length) {
                        printHelp();
                        return null;
                    }
                    p.modelPath = args[++i];
                }
                case "--format" -> {
                    if (i + 1 >= args.length) {
                        printHelp();
                        return null;
                    }
                    p.format = args[++i];
                }
                case "--seq" -> {
                    if (i + 1 >= args.length) {
                        printHelp();
                        return null;
                    }
                    p.fastaPath = args[++i];
                }
                case "--maf" -> {
                    if (i + 1 >= args.length) {
                        printHelp();
                        return null;
                    }
                    p.mafPath = args[++i];
                }
                case "--probabilities" -> p.probabilities = true;
                default -> {

                }
            }
        }

        if (p.modelPath == null || p.format == null) {
            printHelp();
            return null;
        }
        if (!p.format.equals("txt") && !p.format.equals("html")) {
            printHelp();
            return null;
        }

        boolean hasSeq = (p.fastaPath != null);
        boolean hasMaf = (p.mafPath != null);

        //  --seq or --maf
        if (hasSeq == hasMaf) {
            printHelp();
            return null;
        }

        return p;
    }

    private static boolean hasArg(String[] args, String key) {
        for (String a : args) if (a.equals(key)) return true;
        return false;
    }

    private static void printHelp() {
        System.out.println("Training:");
        System.out.println("  java -jar gor4.jar --db <seclib-file> --model <model-file>");
        System.out.println();
        System.out.println("Prediction:");
        System.out.println("  java -jar gor4.jar [--probabilities] --model <model-file> --format <txt|html> {--seq <fasta file>|--maf <multiple-alignment-folder>}");
    }

    static class TrainArgs {
        String dbPath;
        String method;
        String modelPath;
    }

    static class PredictArgs {
        String modelPath;
        String format;          // txt or html
        String fastaPath;       // --seq
        String mafPath;         // --maf
        boolean probabilities;
    }


    private static char probToDigit(double p) {
        if (p < 0) p = 0;
        if (p > 1) p = 1;
        int d = (int) Math.round(p * 9.0);
        if (d < 0) d = 0;
        if (d > 9) d = 9;
        return (char) ('0' + d);
    }


}