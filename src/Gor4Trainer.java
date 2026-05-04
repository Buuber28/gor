import java.util.Map;
import java.util.List;

public class Gor4Trainer {

    private final int halfWindow = 8;
    private final int windowSize = 17;


    // index mapping for each AA
    private final Map<Character, Integer> aaIndex;

    public Gor4Trainer(Map<Character, Integer> aaIndex) {
        this.aaIndex = aaIndex;
    }

    public Gor4Model train(List<SeclibEntry> data){


        Gor4Model model = new Gor4Model(halfWindow, aaIndex);

        //training loop, goes through each protein seq in data
        for (int i = 0; i < data.size(); i++){
            String seq = data.get(i).as;
            String ss = data.get(i).ss;

            // skip bad entries
            if (seq == null || ss == null) continue;
            if (seq.length() != ss.length()) continue;


            //iterate center positions, but only pos with full window
            for(int center = halfWindow;
                center < seq.length() - halfWindow;
                center++){

                char centerAA = seq.charAt(center);
                Integer centerAAIdx = aaIndex.get(Character.toUpperCase(centerAA));
                if (centerAAIdx == null) continue; // unknown aa

                char centerState = ss.charAt(center);
                if (centerState != 'C' && centerState != 'E' && centerState != 'H') continue;

                //getting all the aa pairs using 0-16 indexing
                for (int p = 0; p < windowSize;p++){
                    int posP = center - halfWindow + p;
                    char aaP = seq.charAt(posP);

                    Integer aaPIdx = aaIndex.get(Character.toUpperCase(aaP));
                    if (aaPIdx == null) continue;

                    // 4D  counts
                    model.addCount4D(centerState, centerAAIdx, p, aaPIdx);

                    for (int q = p + 1; q < windowSize; q++ ){
                        int posQ = center - halfWindow +q;
                        char aaQ = seq.charAt(posQ);

                        Integer aaQIdx = aaIndex.get(Character.toUpperCase(aaQ));
                        if (aaQIdx == null) continue;
                        model.addCount6D(centerState, centerAAIdx, p, aaPIdx, q, aaQIdx);
                    }


                }

            }



        }


        return model;

    }



}
