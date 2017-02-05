import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

public class DataAnalysis {
	static double ssr;
	
	public static void main(String[] args) throws IOException { 
//		writeDifferential();
		ssres();
//		logireg();
//		linreg();
//		powreg();
	}
	
	/**
	 * y = a((x-b)/m)^h + c
	 * @throws IOException 
	 */
	static void powreg() throws IOException {
		double a = 0.5;
		double b = 1470;
		double c = 2.4;
		double m = 2212.97;
		double h = 0.45;
		double[] r = new double[7];
		for (int i=-3; i<=3; i++) {
			r[i+3] = powr(a, b+i*0.001, c, m, h);
		}
		b += (findMinIndex(r)-3)*0.001;
		for (int i=-3; i<=3; i++) {
			r[i+3] = powr(a, b, c, m+i*0.001, h);
		}
		m += (findMinIndex(r)-3)*0.001;
		for (int i=-3; i<=3; i++) {
			r[i+3] = powr(a, b, c+i*0.001, m, h);
		}
		c += (findMinIndex(r)-3)*0.001;
		for (int i=-3; i<=3; i++) {
			r[i+3] = powr(a, b, c, m, h+i*0.001);
		}
		h += (findMinIndex(r)-3)*0.001;
		for (int i=-3; i<=3; i++) {
			r[i+3] = powr(a+i*0.001, b, c, m, h);
		}
		a += (findMinIndex(r)-3)*0.001;
		BufferedWriter w = new BufferedWriter(new FileWriter("regression.txt"));
		w.write("Power Regression\n" + "y = " + a + " * ((x - " + b + ")/" + m + ")^" + h + " + " + c + "\n");
		w.write("r^2 = " + r[findMinIndex(r)]);
		w.close();
	}
	
	static double powr(double a, double b, double c, double m, double h) throws FileNotFoundException {
		double result = 1; String[] l; double x; double yi; double yc; double sum = 0;
		Scanner s = new Scanner(new File("diff.txt"));
		while (s.hasNextLine()) {
			l = s.nextLine().split("    ");
			x = Double.parseDouble(l[0]);
			yi = Double.parseDouble(l[1]);
			yc = a*Math.pow((x-b)/m, h) + c;
			sum += (yi - yc) * (yi - yc);
		}
		result -= sum/ssr;
		return result;
	}
	
	static void linreg() throws IOException { 
        int n = 0;
        double[] x = new double[613];
        double[] y = new double[613];

        // first pass: read in data, compute xbar and ybar
        double sumx = 0.0, sumy = 0.0, sumx2 = 0.0; String[] l;
        Scanner s = new Scanner(new File("diff.txt"));
        while(s.hasNextLine()) {
        	l = s.nextLine().split("    ");
            x[n] = Double.parseDouble(l[0]);
            y[n] = Double.parseDouble(l[1]);
            sumx  += x[n];
            sumx2 += x[n] * x[n];
            sumy  += y[n];
            n++;
        }
        double xbar = sumx / n;
        double ybar = sumy / n;

        // second pass: compute summary statistics
        double xxbar = 0.0, yybar = 0.0, xybar = 0.0;
        for (int i = 0; i < n; i++) {
            xxbar += (x[i] - xbar) * (x[i] - xbar);
            yybar += (y[i] - ybar) * (y[i] - ybar);
            xybar += (x[i] - xbar) * (y[i] - ybar);
        }
        double beta1 = xybar / xxbar;
        double beta0 = ybar - beta1 * xbar;
        // analyze results
        int df = n - 2;
        double rss = 0.0;      // residual sum of squares
        double ssr = 0.0;      // regression sum of squares
        for (int i = 0; i < n; i++) {
            double fit = beta1*x[i] + beta0;
            rss += (fit - y[i]) * (fit - y[i]);
            ssr += (fit - ybar) * (fit - ybar);
        }
        double R2    = ssr / yybar;
        double svar  = rss / df;
        double svar1 = svar / xxbar;
        double svar0 = svar/n + xbar*xbar*svar1;
        System.out.println("y   = " + beta1 + " * x + " + beta0);
        System.out.println("R^2                 = " + R2);
        BufferedWriter b = new BufferedWriter(new FileWriter("regression.txt"));
        b.write("Linear Regression\n" + "y  = " + beta1 + " * x + " + beta0 + "\n");
        b.write("r^2 = " + R2);
        b.close(); s.close();
    }
	
	static void logireg() throws IOException {
		double c = 0.242;
		double a = 1717.69;
		double m = 36.05;
		double h = 2.45;
		double[] r = new double[7];
		for (int i=-3; i<=3; i++) {
			r[i+3] = logir(c + i*0.001, a, m, h);
		}
		c += (findMinIndex(r)-3)*0.001;
		for (int i=-3; i<=3; i++) {
			r[i+3] = logir(c, a + i*0.001, m, h);
		}
		a += (findMinIndex(r)-3)*0.001;
		for (int i=-3; i<=3; i++) {
			r[i+3] = logir(c, a, m + i*0.001, h);
		}
		m += (findMinIndex(r)-3)*0.001;
		for (int i=-3; i<=3; i++) {
			r[i+3] = logir(c, a, m, h + 0.001*i);
		}
		h += (findMinIndex(r)-3)*0.001;
		BufferedWriter w = new BufferedWriter(new FileWriter("regression.txt"));
		w.write("Logistic Regression\n" + "y = " + c + "/(1 + e^-((x - " + a + ")/" + m + ")) + " + h + "\n");
		w.write("r^2 = " + r[findMinIndex(r)]);
		w.close();
	}
	
	private static int findMinIndex(double[] a) {
		int min = 0;
		for (int i=1; i<a.length; i++) {
			if (Math.abs(a[i]) < Math.abs(a[min]))
				min = i;
		}
		return min;
	}
	
	static double logir(double c, double a, double m, double h) throws FileNotFoundException {
		double result = 1; String[] l; double x; double yi; double yc; double sum = 0;
		Scanner s = new Scanner(new File("diff.txt"));
		while (s.hasNextLine()) {
			l = s.nextLine().split("    ");
			x = Double.parseDouble(l[0]);
			yi = Double.parseDouble(l[1]);
			yc = c/(1 + Math.pow(Math.E, -(x - a)/m)) + h;
			sum += (yi - yc) * (yi - yc);
		}
		result -= sum/ssr;
		return result;
	}
	
	static double ssres() throws FileNotFoundException {
		double result = 0; double[] r = new double[615]; String[] l; double yhat = 0;
		Scanner s = new Scanner(new File("diff.txt"));
		while (s.hasNextLine()) {
			l = s.nextLine().split("    ");
			r[Integer.parseInt(l[0]) - 1583] = Double.parseDouble(l[1]);
			yhat += r[Integer.parseInt(l[0]) - 1583]/613.0;
		}
		for (int i=0; i<615; i++) {
			if (r[i] != 0) {
				result += (r[i] - yhat) * (r[i] - yhat);
			}
		}
		s.close(); 
		System.out.println(yhat);
		System.out.print(result);
		ssr = result;
		return result;
	}
	
	/**
	 * 1583 - 2197 date range because precision is bad somehow
	 * </br> 
	 * calculates difference between vsop and newcomb
	 * @throws IOException
	 */
	static void writeDifferential() throws IOException {
		Scanner vs = new Scanner(new File("vsopveq.txt"));
		Scanner nw = new Scanner(new File("newcomb_f.txt"));
		BufferedWriter writer = new BufferedWriter(new FileWriter("diff.txt"));
		String v; String n; String dif; double d; int h; int m; int x = 0;
		while (vs.hasNextLine() && nw.hasNextLine()) {
			v = vs.nextLine(); n = nw.nextLine();
			if (!v.substring(0,4).equals(n.substring(0,4))) {
				System.out.println("Error! Years do not line up - " + v.substring(0,4));
				break;
			}
			d = Integer.parseInt(v.substring(5, 7)) - Integer.parseInt(n.substring(8, 10));
			h = Integer.parseInt(v.substring(8, 10)) - Integer.parseInt(n.substring(11, 13));
			m = Integer.parseInt(v.substring(11)) - Integer.parseInt(n.substring(14));
			while ((m < 0 && h > 0) || (h < 0 && d > 0)) {
				if (m < 0 && h > 0) {
					m += 60;
					h --;
				}
				if (h < 0 && d > 0) {
					h += 24;
					d --;
				}
			}
			d += ((double) h * 60.0 + (double) m)/(24.0 * 60.0);
			if (d < 2)
				d++;
			writer.write(v.substring(0,4) + "    " + d + "\n");
			x++;
		}
		//TODO
		vs.close();
		nw.close();
		writer.close();
	}
}
