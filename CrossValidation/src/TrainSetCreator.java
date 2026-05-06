import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

public class TrainSetCreator {

    public static void main(String[] args) {
        String dbPath = null;
        String outPrefix = null;
        Integer kFolds = null;
        String mafSource = null;
        long seed = 42L;

        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "--db":
                    if (i + 1 >= args.length) {
                        printHelp();
                        return;
                    }
                    dbPath = args[++i];
                    break;

                case "--out":
                    if (i + 1 >= args.length) {
                        printHelp();
                        return;
                    }
                    outPrefix = args[++i];
                    break;
                case "--mafSource":
                    if (i + 1 >= args.length) {
                        printHelp();
                        return;
                    }
                    mafSource = args[++i];
                    break;

                case "--cv":
                    if (i + 1 >= args.length) {
                        printHelp();
                        return;
                    }
                    kFolds = Integer.parseInt(args[++i]);
                    break;

                case "--seed":
                    if (i + 1 >= args.length) {
                        printHelp();
                        return;
                    }
                    seed = Long.parseLong(args[++i]);
                    break;

                default:
                    printHelp();
                    return;
            }
        }

        if (dbPath == null || outPrefix == null || kFolds == null) {
            printHelp();
            return;
        }

        if (kFolds < 2) {
            throw new IllegalArgumentException("--cv must be at least 2");
        }

        List<DatasetEntry> entries = DatasetReader.readSeclib(dbPath);

        if (mafSource != null) {
            entries = filterByMaf(entries, mafSource);
            System.out.println("kept " + entries.size() + " entries with matching alignment files");
        }

        List<DatasetEntry> shuffled = new ArrayList<>(entries);
        Collections.shuffle(shuffled, new Random(seed));

        List<List<DatasetEntry>> folds = CrossValidationSplitter.splitKFolds(shuffled, kFolds);

        for (int i = 0; i < folds.size(); i++) {
            List<DatasetEntry> test = folds.get(i);
            List<DatasetEntry> train = new ArrayList<>();

            for (int j = 0; j < folds.size(); j++) {
                if (j == i) continue;
                train.addAll(folds.get(j));
            }

            String foldName = String.valueOf(i + 1);
            FoldWriter.writeSeclib(train, outPrefix + "_fold" + foldName + "_train.db");
            FoldWriter.writeSeclib(test, outPrefix + "_fold" + foldName + "_test.db");

            System.out.println("fold " + foldName + ": train=" + train.size() + " test=" + test.size());
        }
    }

    // HELPERS
    public static void printHelp() {
        System.out.println("Usage: java TrainSetCreator --db <seclib-file> --out <output-prefix> --cv <k> [--seed <long>] [--mafSource <alignment-folder>]");
        System.out.println("--db        input seclib file with > / AS / SS blocks");
        System.out.println("--out       output prefix for generated fold files");
        System.out.println("--cv        number of folds for cross validation");
        System.out.println("--seed      optional RNG seed (default: 42)");
        System.out.println("--mafSource optional folder with per-protein alignment");
    }

    public static List<DatasetEntry> filterByMaf(List<DatasetEntry> entries, String mafSource) {
        List<DatasetEntry> kept = new ArrayList<>();

        java.io.File sourceDir = new java.io.File(mafSource);
        java.io.File[] sourceFiles = sourceDir.listFiles();
        if (sourceFiles == null) {
            throw new RuntimeException("Could not list files in maf source: " + mafSource);
        }

        for (DatasetEntry e : entries) {
            if (findAlignmentFile(sourceFiles, e.id) != null) {
                kept.add(e);
            }
        }

        return kept;
    }

    public static java.io.File findAlignmentFile(java.io.File[] files, String id) {
        String[] exts = {".aln", ".ma", ".maf"};

        // exact match first
        for (String ext : exts) {
            for (java.io.File f : files) {
                if (!f.isFile()) continue;
                if (f.getName().equals(id + ext)) return f;
            }
        }

        // then unique prefix match
        java.io.File hit = null;
        for (java.io.File f : files) {
            if (!f.isFile()) continue;

            String n = f.getName();
            if (!(n.endsWith(".aln") || n.endsWith(".ma") || n.endsWith(".maf"))) continue;

            if (n.startsWith(id)) {
                if (hit != null) {
                    throw new RuntimeException("Multiple alignment files match id: " + id);
                }
                hit = f;
            }
        }

        return hit;
    }

}
