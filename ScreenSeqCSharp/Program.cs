using System.IO;
using System;
using System.Linq;
using System.IO.Compression;
using PureHDF;
using System.Xml;

public class ScreenSeq {
    public static string[] SSBarcodes() {
        return new string[]{
            "GACGCGCGTTGTCAT", 
            "CGTCCTAGGACATAT",
            "CTCCAGCTCGAGCTC",
            "GGACGCAACTTAAGA",
            "GCGAACCGCCGCTG",
            "GTGCAAGAGTTGGCG",
            "GAAGAAGCGTTATTC",
            "TGGTGACAAGTATCT",
            "CTGGTCGAACTTCGT"};
    }

    public static void Main() {
        var inputPath = "C:/Users/jonoj/Repositories/ScreenSeq/Data/3_8_2023_Validation/BarcodesFASTQ/";
        var cellDataPath = "C:/Users/jonoj/Repositories/ScreenSeq/Data/3_8_2023_Validation/filtered_feature_bc_matrix.h5";
        var outputPath = "C:/Users/jonoj/Repositories/ScreenSeq/Data/3_8_2023_Validation/ssCounts_HD2.csv";

        var cellBarcodes = LoadCellBarcodes(cellDataPath);
        var ssBarcodes = SSBarcodes();
        var cellBCCounts = cellBarcodes.ToDictionary(cellBarcode => cellBarcode,
            cellBarcode => ssBarcodes.ToDictionary(ssBarcode => ssBarcode,
                                                    ssBarcode => new HashSet<string>()));
        Dictionary<string, string> cellBcCorrectionMap = new Dictionary<string, string>();
        Dictionary<string, string> ssBcCorrectionMap = new Dictionary<string, string>();

        var r1Files = Directory.EnumerateFiles(inputPath, "*R1*").ToArray();
        var r2Files = Directory.EnumerateFiles(inputPath, "*R2*").ToArray();

        int readNumber = 0;
        int lowQuality = 0;
        int badBarcodes = 0;
        int distCutoff = 2;
        for (int laneIndex = 0; laneIndex < r1Files.Length; laneIndex++) {
            using var compressedStream1 = File.Open(r1Files[laneIndex], FileMode.Open);
            using var compressedStream2 = File.Open(r2Files[laneIndex], FileMode.Open);
            using var decompressor1 = new GZipStream(compressedStream1, CompressionMode.Decompress);
            using var decompressor2 = new GZipStream(compressedStream2, CompressionMode.Decompress);
            using var r1Reader = new StreamReader(decompressor1);
            using var r2Reader = new StreamReader(decompressor2);

            string read1Identifier, read2Identifier, read1Quality, read2Quality, read1, read2 = "";
            read1Identifier = read2Identifier = read1Quality = read2Quality = read1 = read2;

            while (NextRead(r1Reader, r2Reader,
                ref read1Identifier, ref read1Quality, ref read1, 
                ref read2Identifier, ref read2Quality, ref read2)) {
                readNumber++;
                if (readNumber % 100000 == 0) {
                    Console.WriteLine($"Completed read {{{readNumber}}}. {{{lowQuality}}} low-quality reads detected so far. {{{readNumber - badBarcodes}}} good barcodes detected.");
                }

                if (!IdentifiersMatch(read1Identifier, read2Identifier)) {
                    Console.WriteLine($"WARNING: identifiers do not match for read {{{readNumber}}}.");
                }
                if(MinQuality(read1Quality+read2Quality) < 14) {
                    lowQuality++;
                    continue;
                }

                var ssBarcode = read2.Substring(10, 15);
                if(!ssBcCorrectionMap.ContainsKey(ssBarcode)) {
                    var bc = GetMostSimilarBarcode(ssBarcode, ssBarcodes, out int dist);
                    ssBcCorrectionMap[ssBarcode] = dist <= distCutoff ? bc : "";
                }
                ssBarcode = ssBcCorrectionMap[ssBarcode];

                var cellBarcode = read1.Substring(0, 16);
                if (!cellBcCorrectionMap.ContainsKey(cellBarcode)) {
                    var bc = GetMostSimilarBarcode(cellBarcode, cellBarcodes, out int dist);
                    cellBcCorrectionMap[cellBarcode] = dist <= distCutoff ? bc : "";
                }
                cellBarcode = cellBcCorrectionMap[cellBarcode];


                if (ssBarcode == "" || cellBarcode == "") {
                    badBarcodes++;
                    continue;
                }

                var umi = read1.Substring(16, 12);
                cellBCCounts[cellBarcode][ssBarcode].Add(umi);
            }
        }

        List<string> data = new List<string> {
            "CellBC,SS4,SS5,SS6,SS7,SS8,SS9,SS10,SS11,SSN"
        };
        Console.WriteLine("Summing...");
        foreach (var cellBarcode in cellBarcodes) {
            var counts = ssBarcodes.Select(ssBarcode => cellBCCounts[cellBarcode][ssBarcode].Count);
            string newLine = cellBarcode + "," + string.Join(",", counts);
            data.Add(newLine);
        }
        Console.WriteLine("Saving...");
        File.WriteAllLines(outputPath, data);
    }

    public static int HammingDistance(string a, string b) {
        int wrong = 0;
        for(int i = 0; i < a.Length; i++) {
            if (a[i] != b[i]) {
                wrong++;
            }
        }
        return wrong;
    }


    public static string GetMostSimilarBarcode(string inBarcode, string[] barcodes, out int minDistance) {
        minDistance = 1000;
        int mostSimilar = -1;
        for(int i = 0; i < barcodes.Length; i++) {
            var dist = HammingDistance(barcodes[i], inBarcode);
            if(dist < minDistance) {
                minDistance = dist;
                mostSimilar = i;
            }
        }
        return barcodes[mostSimilar];
    }

    public static Dictionary<char, char> invert = new Dictionary<char, char>() { 
        { 'A', 'T' } ,
        { 'T', 'A' },
        {'G', 'C' },
        {'C', 'G' }
    };
    public static string FlipCellBC(string bc) {
        return bc.Substring(0, 7) + invert[bc[7]] + invert[bc[8]] + bc.Substring(9,7);
    }

    public static string[] LoadCellBarcodes(string path) {
        var file = H5File.OpenRead(path);
        var dataset = file.Dataset("/matrix/barcodes").ReadString();
        return dataset.Select(x => FlipCellBC(x.Substring(0, 16))).ToArray();
    }

    public static int MinQuality(string s) {
        return s.Min(x => (int)x) - 33;
    }

    public static bool IdentifiersMatch(string i1, string i2) {
        return i1.Split(" ")[0] == i2.Split(" ")[0];
    }

    public static bool NextRead(StreamReader r1Reader, StreamReader r2Reader, 
        ref string i1, ref string q1, ref string s1,
        ref string i2, ref string q2, ref string s2) {
        for (int i = 0; i < 4; i++) {
            var r1 = r1Reader.ReadLine();
            var r2 = r2Reader.ReadLine();
            if(r1 == null || r2 == null) {
                return false;
            }
            if(i == 0) {
                i1 = r1;
                i2 = r2;
            } else if(i == 1) {
                s1 = r1;
                s2 = r2;
            } else if(i == 3) {
                q1 = r1;
                q2 = r2;
            }
        }
        return true;
    }
}