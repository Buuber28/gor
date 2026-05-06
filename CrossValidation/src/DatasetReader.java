import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class DatasetReader {

    public static List<DatasetEntry> readSeclib(String path) {
        List<DatasetEntry> entries = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String line;
            DatasetEntry current = null;

            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                if (line.startsWith(">")) {
                    if (current != null) {
                        entries.add(current);
                    }
                    current = new DatasetEntry();
                    current.id = line.substring(1).trim();
                } else if (line.startsWith("AS")) {
                    if (current == null) continue;
                    current.as = line.substring(2).trim();
                } else if (line.startsWith("SS")) {
                    if (current == null) continue;
                    current.ss = line.substring(2).trim();
                }
            }

            if (current != null) {
                entries.add(current);
            }



        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        return entries;
    }

}
