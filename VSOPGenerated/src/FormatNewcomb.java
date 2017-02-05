import java.util.Scanner;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;

public class FormatNewcomb {
	public static void main(String[] args) throws IOException {
		format("newcomb_unf.txt");
	}
	
	public static void format(String filename) throws IOException {
		Scanner s = new Scanner(new File(filename));
		String l;
		String w;
		BufferedWriter writer = new BufferedWriter(new FileWriter("newcomb_f.txt"));
		while (s.hasNextLine()) {
			l = s.nextLine();
			if (l.equals("")) {
				l = s.nextLine();
			}
			w = l.substring(0, l.indexOf("    ", l.indexOf("    ") + 1));
			writer.write(w + "\n");
		}
		s.close();
		writer.close();
	}
}
