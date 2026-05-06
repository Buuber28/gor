import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class Validation {

    //prediction file reading
    //expects:
    //>id...
    //AS...
    //PS...
    //will ignore probabilities if present
    public static List<ValidationEntry> readPredictionFile(String path) {
        List<ValidationEntry> entries = new ArrayList<>();

        try (BufferedReader br = openReader(path)) {
            String line;
            ValidationEntry current = null;
            String activeField = null;

            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                if (line.startsWith(">")) {
                    if (current != null) {
                        entries.add(current);
                    }

                    current = new ValidationEntry();
                    current.id = line.substring(1).trim();
                    activeField = null;

                } else if (line.startsWith("AS")) {
                    if (current == null) continue;
                    current.as = line.substring(2).trim();
                    activeField = "AS";

                } else if (line.startsWith("PS")) {
                    if (current == null) continue;
                    current.ps = line.substring(2).trim();
                    activeField = "PS";

                } else {
                    if (current == null || activeField == null) continue;

                    if (activeField.equals("AS")) {
                        current.as += line;
                    } else if (activeField.equals("PS")) {
                        current.ps += line;
                    }
                }

                // ignore PH / PE / PC probability lines and other unsupported lines
            }

            if (current != null) {
                entries.add(current);
            }

        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        return entries;
    }


    //reference file reading
    //expects:
    //<id...
    //AS...
    //SS...
    public static List<ValidationEntry> readReferenceFile(String path) {
        List<ValidationEntry> entries = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(path))) {

            String line;

            ValidationEntry current = null;


            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                if (line.startsWith(">")) {
                    //flush previous
                    if (current != null) {
                        entries.add(current);
                    }

                    current = new ValidationEntry();
                    current.id = line.substring(1).trim();

                } else if (line.startsWith("AS")) {
                    if (current == null) continue;
                    current.as = line.substring(2).trim();

                } else if (line.startsWith("SS")) {
                    if (current == null) continue;
                    ;
                    current.ss = line.substring(2).trim();
                }

                // ignore AS and other lines

            }

            // flush last entry
            if (current != null) {
                entries.add(current);
            }

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return entries;

    }


    //tries to match predictions and references based on the ID
    public static List<ValidationEntry> matchIds(List<ValidationEntry> prediction, List<ValidationEntry> reference) {
        List<ValidationEntry> matchedEntries = new ArrayList<>();

        for (ValidationEntry pred : prediction) {
            for (ValidationEntry ref : reference) {
                if (pred.id != null && pred.id.equals(ref.id)) {
                    ValidationEntry merged = new ValidationEntry();
                    merged.id = pred.id;
                    merged.ss = ref.ss;
                    merged.ps = pred.ps;
                    merged.as = ref.as; //now null because reader ignores AS

                    matchedEntries.add(merged);
                    break;
                }
            }
        }
        return matchedEntries;

    }

    //computes overall q3 and for each state and stores them in an array [q3, qH, qE, qC]
    public static double[] computeQ3(List<ValidationEntry> matchedId) {

        int correct = 0;
        int total = 0;

        int correctH = 0;
        int correctE = 0;
        int correctC = 0;

        int totalH = 0;
        int totalE = 0;
        int totalC = 0;


        for (ValidationEntry e : matchedId) {
            if (e.ps == null || e.ss == null) continue;

            int len = Math.min(e.ps.length(), e.ss.length());

            for (int i = 0; i < len; i++) {
                char pred = e.ps.charAt(i);
                char real = e.ss.charAt(i);

                //ignore edges and gaps
                if (pred == '-') continue;

                total++;

                if (pred == real) correct++;

                if (real == 'H') {
                    totalH++;
                    if (pred == 'H') {
                        correctH++;
                    }
                } else if (real == 'E') {
                    totalE++;
                    if (pred == 'E') {
                        correctE++;
                    }
                } else if (real == 'C') {
                    totalC++;
                    if (pred == 'C') {
                        correctC++;
                    }
                }


            }
        }

        double q3 = total == 0 ? Double.NaN : 100.0 * correct / total;
        double qH = totalH == 0 ? Double.NaN : 100.0 * correctH / totalH;
        double qE = totalE == 0 ? Double.NaN : 100.0 * correctE / totalE;
        double qC = totalC == 0 ? Double.NaN : 100.0 * correctC / totalC;

        return new double[]{q3, qH, qE, qC};
    }

    public static double[] computeSOV(List<ValidationEntry> matchedId) {

        double numAll = 0.0, denAll = 0.0;
        double numH = 0.0, denH = 0.0;
        double numE = 0.0, denE = 0.0;
        double numC = 0.0, denC = 0.0;

        for (ValidationEntry e : matchedId) {
            if (e.ps == null || e.ss == null) continue;
            if (e.ps.isEmpty() || e.ss.isEmpty()) continue;

            String[] filtered = filterForSOV(e.ss, e.ps);
            String ssFiltered = filtered[0];
            String psFiltered = filtered[1];

            if (ssFiltered.isEmpty() || psFiltered.isEmpty()) continue;

            List<Segment> refSegs = extractSegments(ssFiltered);
            List<Segment> predSegs = extractSegments(psFiltered);

            for (Segment ref : refSegs) {
                char state = ref.state;

                if (state != 'H' && state != 'E' && state != 'C') continue;

                double segLen = ref.length();
                boolean hadOverlap = false;

                for (Segment pred : predSegs) {
                    if (pred.state != state) continue;

                    int[] vals = ovMinMaxDelta(ref, pred);
                    int minov = vals[0];
                    int maxov = vals[1];
                    int delta = vals[2];

                    if (minov == 0) continue;

                    hadOverlap = true;
                    double contribution = ((double) (minov + delta) / maxov) * segLen;

                    numAll += contribution;
                    denAll += segLen;

                    switch (state) {
                        case 'H':
                            numH += contribution;
                            denH += segLen;
                            break;
                        case 'E':
                            numE += contribution;
                            denE += segLen;
                            break;
                        case 'C':
                            numC += contribution;
                            denC += segLen;
                            break;
                    }
                }

                // S'(i)  =  observed segment with no overlap
                if (!hadOverlap) {
                    denAll += segLen;

                    switch (state) {
                        case 'H':
                            denH += segLen;
                            break;
                        case 'E':
                            denE += segLen;
                            break;
                        case 'C':
                            denC += segLen;
                            break;
                    }
                }
            }
        }

        double sovAll = denAll == 0 ? Double.NaN : 100.0 * numAll / denAll;
        double sovH   = denH   == 0 ? Double.NaN : 100.0 * numH   / denH;
        double sovE   = denE   == 0 ? Double.NaN : 100.0 * numE   / denE;
        double sovC   = denC   == 0 ? Double.NaN : 100.0 * numC   / denC;

        return new double[]{sovAll, sovH, sovE, sovC};
    }
    public static void main(String[] args) {
        String predPath = null;
        String refPath = null;
        String sumPath = null;
        String detPath = null;
        String outFormat = null;

        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "-p":
                    if (i + 1 >= args.length) {
                        printHelp();
                        return;
                    }
                    predPath = args[++i];
                    break;

                case "-r":
                    if (i + 1 >= args.length) {
                        printHelp();
                        return;
                    }
                    refPath = args[++i];
                    break;

                case "-s":
                    if (i + 1 >= args.length) {
                        printHelp();
                        return;
                    }
                    sumPath = args[++i];
                    break;

                case "-d":
                    if (i + 1 >= args.length) {
                        printHelp();
                        return;
                    }
                    detPath = args[++i];
                    break;

                case "-f":
                    if (i + 1 >= args.length) {
                        printHelp();
                        return;
                    }
                    outFormat = args[++i];
                    break;

                default:
                    printHelp();
                    return;
            }
        }

        if (predPath == null || refPath == null || sumPath == null || detPath == null || outFormat == null) {
            printHelp();
            return;
        }

        if (!outFormat.equals("txt") && !outFormat.equals("html")) {
            printHelp();
            return;
        }

        List<ValidationEntry> prediction = readPredictionFile(predPath);
        List<ValidationEntry> reference = readReferenceFile(refPath);
        List<ValidationEntry> matched = matchIds(prediction, reference);

        double[] q = computeQ3(matched);
        double[] sov = computeSOV(matched);

        if (outFormat.equals("txt")) {
            writeDetailTxt(matched, detPath);
            writeSummaryTxt(matched, sumPath);
        } else if (outFormat.equals("html")) {
            writeDetailHtml(matched, detPath);
            writeSummaryHtml(matched, sumPath);
        }


    }


    //HELPERS

    public static String fmtScore(double x) {
        if (Double.isNaN(x)) return "-";
        return String.format(java.util.Locale.US, "%.1f", x);
    }

    public static void writeSummaryLine(PrintWriter w, String label, List<Double> values) {
        double[] s = fiveNumberSummary(values);

        w.println(
                label + " " +
                        fmtScore(s[0]) + " " +
                        fmtScore(s[1]) + " " +
                        fmtScore(s[2]) + " " +
                        fmtScore(s[3]) + " " +
                        fmtScore(s[4])
        );
    }

    public static void writeSummaryTxt(List<ValidationEntry> matched, String path) {
        try (PrintWriter w = new PrintWriter(path)) {

            List<Double> q3s = new ArrayList<>();
            List<Double> qHs = new ArrayList<>();
            List<Double> qEs = new ArrayList<>();
            List<Double> qCs = new ArrayList<>();

            List<Double> sovs = new ArrayList<>();
            List<Double> sovHs = new ArrayList<>();
            List<Double> sovEs = new ArrayList<>();
            List<Double> sovCs = new ArrayList<>();

            for (ValidationEntry e : matched) {
                List<ValidationEntry> one = new ArrayList<>();
                one.add(e);

                double[] q = computeQ3(one);     // [q3, qH, qE, qC]
                double[] sov = computeSOV(one);  // [sov, sovH, sovE, sovC]


                if (!Double.isNaN(q[0])) q3s.add(q[0]);
                if (!Double.isNaN(q[1])) qHs.add(q[1]);
                if (!Double.isNaN(q[2])) qEs.add(q[2]);
                if (!Double.isNaN(q[3])) qCs.add(q[3]);

                if (!Double.isNaN(sov[0])) sovs.add(sov[0]);
                if (!Double.isNaN(sov[1])) sovHs.add(sov[1]);
                if (!Double.isNaN(sov[2])) sovEs.add(sov[2]);
                if (!Double.isNaN(sov[3])) sovCs.add(sov[3]);
            }


            writeSummaryLine(w, "q3:", q3s);
            writeSummaryLine(w, "qObs_H:", qHs);
            writeSummaryLine(w, "qObs_E:", qEs);
            writeSummaryLine(w, "qObs_C:", qCs);

            writeSummaryLine(w, "SOV:", sovs);
            writeSummaryLine(w, "SOV_H:", sovHs);
            writeSummaryLine(w, "SOV_E:", sovEs);
            writeSummaryLine(w, "SOV_C:", sovCs);

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static void writeDetailTxt(List<ValidationEntry> matched, String path) {

        try (PrintWriter w = new PrintWriter(path)) {

            for (ValidationEntry e : matched) {
                w.println();
                List<ValidationEntry> one = new ArrayList<>();
                one.add(e);

                double[] q = computeQ3(one);
                double[] sov = computeSOV(one);

                w.println(">" + e.id + " "
                        + fmtField(q[0]) + " "
                        + fmtField(sov[0]) + " "
                        + fmtField(q[1]) + " "
                        + fmtField(q[2]) + " "
                        + fmtField(q[3]) + " "
                        + fmtField(sov[1]) + " "
                        + fmtField(sov[2]) + " "
                        + fmtField(sov[3]));



                if (e.as != null) {
                    w.println("AS " + e.as);
                }
                if (e.ps != null) {
                    w.println("PS " + e.ps);
                }
                if (e.ss != null) {
                    w.println("SS " + e.ss);
                }
                w.println();

            }

        } catch (IOException ex) {
            throw new RuntimeException(ex);
        }

    }

    public static void printHelp() {

        System.out.println("Usage: java -jar validateGor.jar -p <predictions> -r <seclib-file> -s <summary file> -d <detailed file> -f <txt|html>"
        );


    }

    public static void writeDetailHtml(List<ValidationEntry> matched, String path) {
        try (PrintWriter w = new PrintWriter(path)) {
            w.println("<html><body><pre>");

            for (ValidationEntry e : matched) {
                w.println();
                List<ValidationEntry> one = new ArrayList<>();
                one.add(e);

                double[] q = computeQ3(one);
                double[] sov = computeSOV(one);

                w.println(">" + e.id + " "
                        + fmtField(q[0]) + " "
                        + fmtField(sov[0]) + " "
                        + fmtField(q[1]) + " "
                        + fmtField(q[2]) + " "
                        + fmtField(q[3]) + " "
                        + fmtField(sov[1]) + " "
                        + fmtField(sov[2]) + " "
                        + fmtField(sov[3]));


                if (e.as != null) w.println("AS " + e.as);
                if (e.ps != null) w.println("PS " + e.ps);
                if (e.ss != null) w.println("SS " + e.ss);
                w.println();
            }

            w.println("</pre></body></html>");
        } catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    public static void writeSummaryHtml(List<ValidationEntry> matched, String path) {
        try (PrintWriter w = new PrintWriter(path)) {

            List<Double> q3s = new ArrayList<>();
            List<Double> qHs = new ArrayList<>();
            List<Double> qEs = new ArrayList<>();
            List<Double> qCs = new ArrayList<>();

            List<Double> sovs = new ArrayList<>();
            List<Double> sovHs = new ArrayList<>();
            List<Double> sovEs = new ArrayList<>();
            List<Double> sovCs = new ArrayList<>();

            for (ValidationEntry e : matched) {
                List<ValidationEntry> one = new ArrayList<>();
                one.add(e);

                double[] q = computeQ3(one);
                double[] sov = computeSOV(one);

                if (!Double.isNaN(q[0])) q3s.add(q[0]);
                if (!Double.isNaN(q[1])) qHs.add(q[1]);
                if (!Double.isNaN(q[2])) qEs.add(q[2]);
                if (!Double.isNaN(q[3])) qCs.add(q[3]);

                if (!Double.isNaN(sov[0])) sovs.add(sov[0]);
                if (!Double.isNaN(sov[1])) sovHs.add(sov[1]);
                if (!Double.isNaN(sov[2])) sovEs.add(sov[2]);
                if (!Double.isNaN(sov[3])) sovCs.add(sov[3]);
            }

            w.println("<html><body><pre>");
            writeSummaryLine(w, "q3:", q3s);
            writeSummaryLine(w, "qObs_H:", qHs);
            writeSummaryLine(w, "qObs_E:", qEs);
            writeSummaryLine(w, "qObs_C:", qCs);

            writeSummaryLine(w, "SOV:", sovs);
            writeSummaryLine(w, "SOV_H:", sovHs);
            writeSummaryLine(w, "SOV_E:", sovEs);
            writeSummaryLine(w, "SOV_C:", sovCs);
            w.println("</pre></body></html>");

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    //returns an array with min overlap, maximum overlap and delta in it,

    public static double[] fiveNumberSummary(List<Double> values) {
        if (values.isEmpty()) return new double[]{0, 0, 0, 0, 0};

        List<Double> sorted = new ArrayList<>(values);
        sorted.sort(Double::compareTo);

        int n = sorted.size();

        double min = sorted.get(0);
        double q25 = sorted.get((int)Math.floor(0.25 * (n - 1)));
        double median = sorted.get((int)Math.floor(0.50 * (n - 1)));
        double q75 = sorted.get((int)Math.floor(0.75 * (n - 1)));
        double max = sorted.get(n - 1);

        return new double[]{min, q25, median, q75, max};
    }


    //returns an array with min overlap, maximum overlap and delta in it,

    public static int[] ovMinMaxDelta(Segment s1, Segment s2) {
        int startmax = Math.max(s1.start, s2.start);
        int endmin = Math.min(s1.end, s2.end);
        int[] result = new int[3];
        if (endmin < startmax) {
            result[0] = 0;// no overlap

        } else {
            result[0] = endmin - startmax + 1;
        }


        int startmin = Math.min(s1.start, s2.start);
        int endmax = Math.max(s1.end, s2.end);

        result[1] = endmax - startmin + 1;
        int d1 = result[1] - result[0];
        int d2 = result[0];
        int d3 = s1.length() / 2;
        int d4 = s2.length() / 2;

        result[2] = Math.min(Math.min(d1, d2), Math.min(d3, d4));

        return result;
    }

    public static List<Segment> extractSegments(String s) {
        List<Segment> segments = new ArrayList<>();

        if (s == null || s.isEmpty()) return segments;

        char currentState = s.charAt(0);
        int start = 0;

        for (int i = 1; i < s.length(); i++) {
            char c = s.charAt(i);


            //new segment
            if (c != currentState) {
                //ignore "-" segments
                if (currentState != '-') {
                    segments.add(new Segment(currentState, start, i - 1));
                }

                currentState = c;
                start = i;
            }
        }

        //flush last
        if (currentState != '-') {
            segments.add(new Segment(currentState, start, s.length() - 1));
        }
        return segments;
    }

    public static class ValidationEntry {
        String id;
        String as;
        String ss;
        String ps;

    }

    public static class Segment {
        char state;
        int start;
        int end;

        public Segment(char state, int start, int end) {
            this.state = state;
            this.start = start;
            this.end = end;
        }

        public int length() {
            return end - start + 1;
        }
    }

    public static BufferedReader openReader(String path) throws IOException {
        if (path.equals("-")) {
            return new BufferedReader(new InputStreamReader(System.in));
        }
        return new BufferedReader(new FileReader(path));
    }

    public static String fmtField(double x) {
        return String.format(java.util.Locale.US, "%5s", fmtScore(x));
    }

    public static String[] filterForSOV(String ss, String ps) {
        StringBuilder ssFiltered = new StringBuilder();
        StringBuilder psFiltered = new StringBuilder();

        int len = Math.min(ss.length(), ps.length());

        for (int i = 0; i < len; i++) {
            char pred = ps.charAt(i);
            char real = ss.charAt(i);

            if (pred == '-') continue; // ignore positions that were not predicted

            ssFiltered.append(real);
            psFiltered.append(pred);
        }

        return new String[]{ssFiltered.toString(), psFiltered.toString()};
    }

}
