import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

public class FoldWriter {

    public static void writeSeclib(List<DatasetEntry> entries, String path) {
        try (PrintWriter w = new PrintWriter(new FileWriter(path))) {
            for (DatasetEntry e : entries) {
                w.println(">" + e.id);
                if (e.as != null) w.println("AS " + e.as);
                if (e.ss != null) w.println("SS " + e.ss);
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
