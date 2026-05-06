import java.util.HashMap;
import java.util.Map;

public class Gor4Model {

    private final int halfWindow;
    private final int windowSize;
    private final Map<Character, Integer> aaIndex;


    private static final double ALPHA = 1.0;//usin a constant to avoid log(0)

    // counts[state][centerAA][aaP][pIdx][aaQ][qIdx]
    private final long[][][][][][] counts6D;

    // counts4D[state][centerAA][aaP][pIdx]
    private final long[][][][] counts4D;

    // weights[state][centerAA][aaP][pIdx][aaQ][qIdx]
    private final double[][][][][][] weights;

    // weights4D[state][centerAA][aaP][pIdx]
    private final double[][][][] weights4D;

    public Gor4Model(int halfWindow, Map<Character, Integer> aaIndex) {
        this.halfWindow = halfWindow;
        this.windowSize = 2 * halfWindow + 1;
        this.aaIndex = aaIndex;

        int states = 3;
        int aa = 20;

        counts6D  = new long[states][aa][aa][windowSize][aa][windowSize];
        counts4D  = new long[states][aa][aa][windowSize];
        weights = new double[states][aa][aa][windowSize][aa][windowSize];
        weights4D = new double[states][aa][aa][windowSize];
    }

    private int stateIndex(char s) {
        if (s == 'C') return 0;
        if (s == 'E') return 1;
        if (s == 'H') return 2;
        return -1;
    }

    //these function is called by the trainer
    public void addCount6D(char centerState, int centerAA, int pIdx, int aaP, int qIdx, int aaQ) {
        int si = stateIndex(centerState);
        if (si < 0) return;
        counts6D[si][centerAA][aaP][pIdx][aaQ][qIdx] += 1L;
    }


    public void addCount4D(char centerState, int centerAA, int pIdx, int aaP) {
        int si = stateIndex(centerState);
        if (si < 0) return;
        counts4D[si][centerAA][aaP][pIdx] += 1L;
    }

    // Copy raw counts from a model file into this model.
    public void loadCounts6D(long[][][][][][] src) {
        for (int s = 0; s < 3; s++) {
            for (int centerAA = 0; centerAA < 20; centerAA++) {
                for (int aaP = 0; aaP < 20; aaP++) {
                    for (int pIdx = 0; pIdx < windowSize; pIdx++) {
                        for (int aaQ = 0; aaQ < 20; aaQ++) {
                            for (int qIdx = 0; qIdx < windowSize; qIdx++) {
                                counts6D[s][centerAA][aaP][pIdx][aaQ][qIdx] = src[s][centerAA][aaP][pIdx][aaQ][qIdx];
                            }
                        }
                    }
                }
            }
        }
    }

    public void loadCounts4D(long[][][][] src) {
        for (int s = 0; s < 3; s++) {
            for (int centerAA = 0; centerAA < 20; centerAA++) {
                for (int aaP = 0; aaP < 20; aaP++) {
                    for (int pIdx = 0; pIdx < windowSize; pIdx++) {
                        counts4D[s][centerAA][aaP][pIdx] = src[s][centerAA][aaP][pIdx];
                    }
                }
            }
        }
    }


    // turns raw counts into weights
    public void finalizeModel() {

        // 4D weights
        // weights4D[state][centerAA][aaP][pIdx]
        for (int centerAA = 0; centerAA < 20; centerAA++) {
            for (int aaP = 0; aaP < 20; aaP++) {
                for (int pIdx = 0; pIdx < windowSize; pIdx++) {


                    for (int s = 0; s < 3; s++) {

                        long countS = counts4D[s][centerAA][aaP][pIdx];

                        long countNotS =
                                counts4D[(s+1)%3][centerAA][aaP][pIdx] +
                                        counts4D[(s+2)%3][centerAA][aaP][pIdx];

                        double num = countS + ALPHA;
                        double den = countNotS + ALPHA;

                        weights4D[s][centerAA][aaP][pIdx] = Math.log(num / den);
                    }
                }
            }
        }

        //weights 6D
        // weights[state][centerAA][aaP][pIdx][aaQ][qIdx]
        for (int centerAA = 0; centerAA < 20; centerAA++) {
            for (int aaP = 0; aaP < 20; aaP++) {
                for (int pIdx = 0; pIdx < windowSize; pIdx++) {
                    for (int aaQ = 0; aaQ < 20; aaQ++) {
                        for (int qIdx = 0; qIdx < windowSize; qIdx++) {



                            for (int s = 0; s < 3; s++) {

                                long countS = counts6D[s][centerAA][aaP][pIdx][aaQ][qIdx];

                                long countNotS =
                                        counts6D[(s+1)%3][centerAA][aaP][pIdx][aaQ][qIdx] +
                                                counts6D[(s+2)%3][centerAA][aaP][pIdx][aaQ][qIdx];

                                double num = countS + ALPHA;
                                double den = countNotS + ALPHA;

                                weights[s][centerAA][aaP][pIdx][aaQ][qIdx] = Math.log(num / den);
                            }

                        }
                    }
                }
            }
        }
    }

    // Score the 3 states (C/E/H) for center position j in seq.
    // Returns a double[3] with a score for each state
    //   score(S) = (2/(2m+1)) * sum_{k,l>k} log odds from Weights6D
    //            - ((2m-1)/(2m+1)) * sum_{k} log odds from Weights4D

    public double[] scoreAt(String seq, int j) {
        double[] scores = new double[3];

        // only use a full window around j
        if (j < halfWindow || j >= seq.length() - halfWindow) {
            return scores; // caller should output '-' at borders anyway
        }

        // center amino acid index
        char centerChar = seq.charAt(j);
        Integer centerAA = aaIndex.get(centerChar);
        if (centerAA == null) {
            // Unknown center AA -> treat as Coil.
            return new double[]{0.0, -1e9, -1e9};
        }



        double pairCoef   = 2.0 / windowSize;              // 2/(2m+1)
        double singleCoef = (windowSize - 2.0) / windowSize; // (2m-1)/(2m+1)

        for (int s = 0; s < 3; s++) {

            double sumSingles = 0.0;
            double sumPairs   = 0.0;

            //  sum over k (4D)
            for (int pIdx = 0; pIdx < windowSize; pIdx++) {
                int posP = j + (pIdx - halfWindow);
                char aaPChar = seq.charAt(posP);
                Integer aaP = aaIndex.get(aaPChar);
                if (aaP == null) continue;

                sumSingles += weights4D[s][centerAA][aaP][pIdx];
            }

            //   sum over k,l>k (6D)
            for (int pIdx = 0; pIdx < windowSize; pIdx++) {
                int posP = j + (pIdx - halfWindow);
                Integer aaP = aaIndex.get(seq.charAt(posP));
                if (aaP == null) continue;

                for (int qIdx = pIdx + 1; qIdx < windowSize; qIdx++) {
                    int posQ = j + (qIdx - halfWindow);
                    Integer aaQ = aaIndex.get(seq.charAt(posQ));
                    if (aaQ == null) continue;

                    sumPairs += weights[s][centerAA][aaP][pIdx][aaQ][qIdx];
                }
            }

            // combine
            scores[s] = pairCoef * sumPairs - singleCoef * sumSingles;
        }

        return scores;
    }



    // Predicts the secondary structure string
    public String predictPS(String seq) {
        StringBuilder ps = new StringBuilder(seq.length());

        // Border positions  are '-'
        for (int j = 0; j < seq.length(); j++) {
            if (j < halfWindow || j >= seq.length() - halfWindow) {
                ps.append('-');
                continue;
            }


            char c = seq.charAt(j);
            if (!aaIndex.containsKey(c)) {
                ps.append("C");   // unknown AA -> predict Coil
                continue;
            }

            double[] sc = scoreAt(seq, j);
            int best = maxIdxDouble(sc);


            String state = switch(best){
                case 0 -> "C";
                case 1 -> "E";
                case 2 -> "H";
                default -> "X"; //error
            };
            ps.append(state);
        }

        return ps.toString();
    }






    //HELPERS

    public int maxIdxDouble(double[] arr){
        int biggest = 0;
        for(int i = 0; i< arr.length;i++){
            if(arr[biggest] < arr[i] ){
                biggest = i;
            }
        }
        return biggest;
    }

    public long[][][][][][] getCounts6D() {
        return counts6D;
    }

    public long[][][][] getCounts4D() {
        return counts4D;
    }

    public double[][][][] getWeights4D() {
        return weights4D;
    }

    public double[][][][][][] getWeights6D() {
        return weights;
    }
}