import java.net.DatagramSocket;
import java.util.ArrayList;
import java.util.List;

public class CrossValidationSplitter {
    public static List<List<DatasetEntry>> splitKFolds(List<DatasetEntry> entries, int k) {
        List<List<DatasetEntry>> folds = new ArrayList<>();
        for (int i = 0; i < k; i++) {
            folds.add(new ArrayList<>());
        }

        for (int i = 0; i < entries.size(); i++) {
            folds.get(i % k).add(entries.get(i));
        }

        return folds;
    }


}
