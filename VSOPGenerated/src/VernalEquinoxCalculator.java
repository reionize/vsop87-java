import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;

public class VernalEquinoxCalculator {

	static final double k = 1.49597870691 / 299792458 / 86400 * Math.pow(10, 11); // days for light to travel 1 AU
	static double hX;
	static double hY;
	static double hZ;
	static double r;
	static double start = 2440587;
	
	
	
	public static void main (String[] args) throws IOException, ParseException {
		double gy = 1452;
		double t;
		BufferedWriter writer = new BufferedWriter(new FileWriter("vsop_unf_1.txt"));
		
//		for (int i=0; i<365; i++) {
//			writer.write(Double.toString(VSOP87C.Earth_X(t)));
//			writer.write(" " + Double.toString((VSOP87C.Earth_Y(t))));
//			t += 1/365350;
//		}
		for (int yy=2300; yy<2320; yy++) {
			String l = yy + "\t" + MarchEquinox(yy);
			if (yy % 5 == 0) {
				System.out.println(yy);
			}
			writer.write(l + "\n");
		}
		writer.close();
	}
	
	
	public static String MarchEquinox (int gy) throws ParseException { 
		
		int[] ymdh = {gy, 3, 20, 12};
		double jcenter = toJulian(ymdh);

	   	double[] a = new double[14];
	   	int i = 0;
	   	double t;

	   	for(double j = jcenter-3; j <= jcenter+3; j++) {
	   		t = (j-2451545)/36525; 
	   		a[i] = geocentric(t, "dec");
	   		a[i+1] = j;
	   		i += 2;
	   	}
	   	
	   	/** $JDEvent  = LaGrange_Interpolate ($DataTable, 0);
	   	$Ymd_HMSTT = JD_To_Ymd_HMS($JDEvent);
	   	$Ymd_HMSTT = UDateTime_To_TZDateTime ($Ymd_HMSTT, $TZhm, 1); */
	   	double jdEvent = LaGrange_Interpolate(a, 0);
	   	String jdF = Double.toString(jdEvent);
	   	
	   	double time = jdEvent - Math.floor(jdEvent);
	   	String hh = Integer.toString((int) Math.floor(time*24));
	   	if (hh.length() < 2) {
	   		hh = "0" + hh;
	   	}
	   	String mm = Integer.toString((int) Math.floor((time*24 - Integer.parseInt(hh))*60));
	   	if (mm.length() < 2) {
	   		mm = "0" + mm;
	   	}
	   	
	   	//start is 01/01/1970
	   	double day = jdEvent - start;
	   	String j = Double.toString(day);
	   	Date date = new SimpleDateFormat("D").parse(j);
		String g = new SimpleDateFormat("dd").format(date);
	   	
	   	return g + " " + hh + " " + mm;
	}
	
	
	static double geocentric (double t, String s) { //use "ra" or "dec"
		hX = VSOP87C.Earth_X(t);
		hY = VSOP87C.Earth_Y(t);
		hZ = VSOP87C.Earth_Z(t);
		r = Math.sqrt(hX * hX + hY * hY + hZ * hZ);
		double ltc = r * k; //this is in days!!
		double ltp = 0;
		double ta = t;
		
		//light time correction
		while (ltc - ltp > 0.000000000001) {
			ltp = ltc;
			ta -= ltc/365250;
			hX = VSOP87C.Earth_X(ta);
			hY = VSOP87C.Earth_Y(ta);
			hZ = VSOP87C.Earth_Z(ta);
			r = Math.sqrt(hX * hX + hY * hY + hZ * hZ);
			ltc = r*k;
		}
		
		double eclLngDeg = Math.toDegrees(Math.atan(hY/hX));
		if (eclLngDeg < 0)
			eclLngDeg += 360;
		double eclLatDeg = Math.toDegrees(Math.atan(-hZ/Math.sqrt(hX * hX + hY * hY)));
		
		//FK5 correction
		double u = Math.toRadians(eclLatDeg);
		double lprime = Math.toRadians(eclLngDeg - 1.297*t- 0.00031*t*t);
		double w = eclLngDeg + (0.03916*(Math.cos(lprime) + Math.sin(lprime))*Math.tan(u) - 0.09033) / 3600;
		
		// Correction for nutation in longitude using the IAU 2000B nutation series
		   eclLngDeg += Delta_Psi_2000B(t);
		   
		// Compute true obliquity of the ecliptic using the IAU 2000B nutation series.
		   double epsTrue = Mean_Epsilon(t) + Delta_Epsilon_2000B(t);
		   double eps = Math.toRadians(epsTrue);
		   u = Math.toRadians(eclLngDeg);
		   w = Math.toRadians(eclLatDeg);
		
		if (s.equals("ra")) {
			double y = Math.sin(u)*Math.cos(eps) - Math.tan(w)*Math.sin(eps);
			double x = Math.cos(u);
			double RADeg = Math.toDegrees(Math.atan(y/x));   if (RADeg < 0) {RADeg += 360;}
			double RAHrs = RADeg / 15;
			return RAHrs;
		}
		
		if (s.equals("dec")) { // for March equinox, dec == 0 is used
			w = Math.asin(Math.sin(w)*Math.cos(eps) + Math.cos(w)*Math.sin(u)*Math.sin(eps));
			double declDeg = Math.toDegrees(w);
			return declDeg;
		}
		return 0.0;
	}
	
	
	public static double toJulian(int[] ymdh) {
		int JGREG= 15 + 31*(10+12*1582);
		int year=ymdh[0];
		int month=ymdh[1]; // jan=1, feb=2,...
		int day=ymdh[2];
		double h=ymdh[3];
		int julianYear = year;
		if (year < 0) julianYear++;
		int julianMonth = month;
		if (month > 2) {
			julianMonth++;
		}
		else {
			julianYear--;
			julianMonth += 13;
		}
				
		double julian = (java.lang.Math.floor(365.25 * julianYear) + java.lang.Math.floor(30.6001*julianMonth) + day + 1720995.0);
		if (day + 31 * (month + 12 * year) >= JGREG) {
			// change over to Gregorian calendar
			int ja = (int)(0.01 * julianYear);
			julian += 2 - ja + (0.25 * ja);
		}
		return java.lang.Math.floor(julian)+h/24;
		}
	
	
	private static double LaGrange_Interpolate (double[] a, double xarg) {

	// Count the total number of data elements. 

	   int totalDataCount = a.length;

	// Number of data pairs.  

	   int n = totalDataCount / 2;
	   
	   double[] X = new double[n];
	   double[] Y = new double[n];

	// Construct separate XY arrays from the array data.

	   for(int i=0; i < totalDataCount; i+=2) {
	       X[i/2] = a[i];
	       Y[i/2] = a[i+1];
//		   System.out.println(a[i]);
	   }

	// Initialize Y summation accumulator.
	   double y = 0.0;
	   
	   double x = xarg;

	// Compute Lagrangian product (Li) for given X argument.

	   for (int i=0; i < n; i++) {
	        double Li = 1.0;

	        for (int j=0; j < n; j++)
	            {
	             if (j != i) // Skip this cycle when j == i
	                {
	                 Li = (Li * (x - X[j])) / (X[i] - X[j]);
	                 
	                }
	   }

//	      Accumulate sum of Yi polynomial terms.

	        y += (Y[i] * Li);

	       }

	   return y;

	}
	
	
	private static double Delta_Epsilon_2000B(double t) {

	   double t2 = t * t;
	   double t3 = t * t2;
	   double t4 = t * t3;

	// -------------------------------------------
	// Compute mean anomaly of the Moon in radians

	   double l  = Math.toRadians((485868.249036 + 1717915923.2178*t + 31.8792*t2
	       + 0.051635*t3 - 0.0002447*t4) / 3600);


	// ------------------------------------------
	// Compute mean anomaly of the Sun in radians

	   double lp = Math.toRadians((1287104.79305 + 129596581.0481*t
	       - 0.5532*t2  + 0.000136*t3 - 0.00001149*t4) / 3600);


	// ------------------------------------------------------------
	// Compute mean argument of the latitude of the Moon in radians

	   double f  = Math.toRadians((335779.526232 + 1739527262.8478*t
	       - 12.7512*t2 - 0.001037*t3 + 0.00000417*t4) / 3600);


	// -----------------------------------------------------------
	// Compute mean elongation of the Moon from the Sun in radians

	   double d  = Math.toRadians((1072260.70369 + 1602961601.2090*t
	       - 6.3706*t2  + 0.006593*t3 - 0.00003169*t4) / 3600);


	// ------------------------------------------------------------------------
	// Compute mean longitude of the mean ascending node of the Moon in radians

	   double om = Math.toRadians((450160.398036 - 6962890.5431*t
	       + 7.4722*t2  + 0.007702*t3 - 0.00005939*t4) / 3600);

	// --------------------------------------------------------------
	// Sum series for nutation in obliquity (dEps) in arc sec * 10E+7

	double s = 0;
	s += (92052331 + 9086*t)*Math.cos(om) + 15377*Math.sin(om);
	s += (5730336 - 3015*t)*Math.cos(2*(f - d + om)) - 4587*Math.sin(2*(f - d + om));
	s += (978459 - 485*t)*Math.cos(2*(f + om)) + 1374*Math.sin(2*(f + om));
	s += (-897492 + 470*t)*Math.cos(2*om) - 291*Math.sin(2*om);
	s += (73871 - 184*t)*Math.cos(lp) - 1924*Math.sin(lp);
	s += (224386 - 677*t)*Math.cos(lp + 2*(f - d + om)) - 174*Math.sin(lp + 2*(f - d + om));
	s -= 6750*Math.cos(l) - 358*Math.sin(l);
	s += (200728 + 18*t)*Math.cos(2*f + om) + 318*Math.sin(2*f + om);
	s += (129025 - 63*t)*Math.cos(l + 2*(f + om)) + 367*Math.sin(l + 2*(f + om));
	s += (-95929 + 299*t)*Math.cos(2*(f - d + om) - lp) + 132*Math.sin(2*(f - d + om) - lp);
	s += (-68982 - 9*t)*Math.cos(2*(f - d) + om) + 39*Math.sin(2*(f - d) + om);
	s += (-53311 + 32*t)*Math.cos(2*(f + om) - l) - 4*Math.sin(2*(f + om) - l);
	s -= 1235*Math.cos(2*d - l) - 82*Math.sin(2*d - l);
	s -= 33228*Math.cos(l + om) + 9*Math.sin(l + om);
	s += 31429*Math.cos(om - l) - 75*Math.sin(om - l);
	s += (25543 - 11*t)*Math.cos(2*(f + d + om) - l) + 66*Math.sin(2*(f + d + om) - l);
	s += 26366*Math.cos(l + 2*f + om) + 78*Math.sin(l + 2*f + om);
	s += (-24236 - 10*t)*Math.cos(2*(f - l) + om) + 20*Math.sin(2*(f - l) + om);
	s -= 1220*Math.cos(2*d) - 29*Math.sin(2*d);
	s += (16452 - 11*t)*Math.cos(2*(f + d + om)) + 68*Math.sin(2*(f + d + om));
	s -= 13870*Math.cos(-2*(lp - f + d - om));
	s += 477*Math.cos(2*(d - l)) - 25*Math.sin(2*(d - l));
	s += (13238 - 11*t)*Math.cos(2*(l + f + om)) + 59*Math.sin(2*(l + f + om));
	s += (-12338 + 10*t)*Math.cos(l + 2*(f - d + om)) - 3*Math.sin(l + 2*(f - d + om));
	s -= 10758*Math.cos(2*f + om - l) + 3*Math.sin(2*f + om - l);
	s -= 609*Math.cos(2*l) - 13*Math.sin(2*l);
	s -= 550*Math.cos(2*f) - 11*Math.sin(2*f);
	s += (8551 - 2*t)*Math.cos(lp + om) - 45*Math.sin(lp + om);
	s -= 8001*Math.cos(-l + 2*d + om) + Math.sin(-l + 2*d + om);
	s += (6850 - 42*t)*Math.cos(2*(lp + f - d + om)) - 5*Math.sin(2*(lp + f - d + om));
	s -= 167*Math.cos(2*(d - f)) - 13*Math.sin(2*(d - f));
	s += 6953*Math.cos(l - 2*d + om) - 14*Math.sin(l - 2*d + om);
	s += 6415*Math.cos(om - lp) + 26*Math.sin(om - lp);
	s += 5222*Math.cos(2*(f + d) + om - l) + 15*Math.sin(2*(f + d) + om - l);
	s += (168 - t)*Math.cos(2*lp) + 10*Math.sin(2*lp);
	s += 3268*Math.cos(l + 2*(f + d + om)) + 19*Math.sin(l + 2*(f + d + om));
	s += 104*Math.cos(2*(f - l)) + 2*Math.sin(2*(f - l));
	s -= 3250*Math.cos(lp + 2*(f + om)) + 5*Math.sin(lp + 2*(f + om));
	s += 3353*Math.cos(2*(f + d) + om) + 14*Math.sin(2*(f + d) + om);
	s += 3070*Math.cos(2*(f + om) - lp) + 4*Math.sin(2*(f + om) - lp);
	s += 3272*Math.cos(2*d + om) + 4*Math.sin(2*d + om);
	s -= 3045*Math.cos(l + 2*(f - d) + om) + Math.sin(l + 2*(f - d) + om);
	s -= 2768*Math.cos(2*(l + f - d + om)) + 4*Math.sin(2*(l + f - d + om));
	s += 3041*Math.cos(2*(d - l) + om) - 5*Math.sin(2*(d - l) + om);
	s += 2695*Math.cos(2*(l + f) + om) + 12*Math.sin(2*(l + f) + om);
	s += 2719*Math.cos(2*(f - d) + om - lp) - 3*Math.sin(2*(f - d) + om - lp);
	s += 2720*Math.cos(om - 2*d) - 9*Math.sin(om - 2*d);
	s -= 51*Math.cos(-l - lp + 2*d) - 4*Math.sin(-l - lp + 2*d);
	s -= 2206*Math.cos(2*(l - d) + om) - Math.sin(2*(l - d) + om);
	s -= 199*Math.cos(l + 2*d) - 2*Math.sin(l + 2*d);
	s -= 1900*Math.cos(lp + 2*(f - d) + om) - Math.sin(lp + 2*(f - d) + om);
	s -= 41*Math.cos(l - lp) - 3*Math.sin(l - lp);
	s += 1313*Math.cos(-2*(l - f - om)) - Math.sin(-2*(l - f - om));
	s += 1233*Math.cos(3*l + 2*(f + om)) + 7*Math.sin(3*l + 2*(f + om));
	s -= 81*Math.cos(-lp + 2*d) - 2*Math.sin(-lp + 2*d);
	s += 1232*Math.cos(l - lp + 2*(f + om)) + 4*Math.sin(l - lp + 2*(f + om));
	s -= 20*Math.cos(d) + 2*Math.sin(d);
	s += 1207*Math.cos(2*(f + d + om) - l - lp) + 3*Math.sin(2*(f + d + om) - l - lp);
	s += 40*Math.cos(2*f - l) - 2*Math.sin(2*f - l);
	s += 1129*Math.cos(-lp + 2*(f + d + om)) + 5*Math.sin(-lp + 2*(f + d + om));
	s += 1266*Math.cos(om - 2*l) - 4*Math.sin(om - 2*l);
	s -= 1062*Math.cos(l + lp + 2*(f + om)) + 3*Math.sin(l + lp + 2*(f + om));
	s -= 1129*Math.cos(2*l + om) + 2*Math.sin(2*l + om);
	s -= 9*Math.cos(lp + d - l);
	s += 35*Math.cos(l + lp) - 2*Math.sin(l + lp);
	s -= 107*Math.cos(l + 2*f) - Math.sin(l + 2*f);
	s += 1073*Math.cos(2*(f - d) + om - l) - 2*Math.sin(2*(f - d) + om - l);
	s += 854*Math.cos(l + 2*om);
	s -= 553*Math.cos(d - l) + 139*Math.sin(d - l);
	s -= 710*Math.cos(2*(f + om) + d) + 2*Math.sin(2*(f + om) + d);
	s += 647*Math.cos(2*(f + 2*d + om) - l) + 4*Math.sin(2*(f + 2*d + om) - l);
	s -= 700*Math.cos(lp + d + om - l);
	s += 672*Math.cos(-2*(lp - f + d) + om);
	s += 663*Math.cos(l + 2*(f + d) + om) + 4*Math.sin(l + 2*(f + d) + om);
	s -= 594*Math.cos(-2*(l - f - d - om)) + 2*Math.sin(-2*(l - f - d - om));
	s -= 610*Math.cos(2*om - l) - 2*Math.sin(2*om - l);
	s -= 556*Math.cos(l + lp + 2*(f - d + om));

	// ----------------------------
	// Return nutation in obliquity expressed in degrees.
	   double rr = 36000000000.0;
	   return s / rr;
	}
	
	
	private static double Mean_Epsilon(double t) {

		t /= 100;

	// ------------------------------------------
	// Compute mean obliquity (e) of the ecliptic expressed in arc seconds.

	   double p = t*t;
	   double e = 84381.448 - 4680.93*t;
	   e -=    1.55*p;  p *= t;
	   e += 1999.25*p;  p *= t;
	   e -=   51.38*p;  p *= t;
	   e -=  249.67*p;  p *= t;
	   e -=   39.05*p;  p *= t;
	   e +=    7.12*p;  p *= t;
	   e +=   27.87*p;  p *= t;
	   e +=    5.79*p;  p *= t;
	   e +=    2.45*p;

	// ---------------------------------------
	// Convert mean obliquity from arc seconds into degrees and return the result.

	   e /= 3600;
	   return e;
	}
	
	
	private static double Delta_Psi_2000B(double t) {

	   double t2 = t*t;
	   double t3 = t*t2;
	   double t4 = t*t3;

	// -------------------------------------------
	// Compute mean anomaly of the Moon in radians

	   double l  = Math.toRadians((485868.249036 + 1717915923.2178*t + 31.8792*t2
	       + 0.051635*t3 - 0.0002447*t4) / 3600);

	// ------------------------------------------
	// Compute mean anomaly of the Sun in radians

	   double lp = Math.toRadians((1287104.79305 + 129596581.0481*t
	       - 0.5532*t2  + 0.000136*t3 - 0.00001149*t4) / 3600);

	// ------------------------------------------------------------
	// Compute mean argument of the latitude of the Moon in radians

	   double f  = Math.toRadians((335779.526232 + 1739527262.8478*t
	       - 12.7512*t2 - 0.001037*t3 + 0.00000417*t4) / 3600);

	// -----------------------------------------------------------
	// Compute mean elongation of the Moon from the Sun in radians

	   double d  = Math.toRadians((1072260.70369 + 1602961601.2090*t
	       - 6.3706*t2  + 0.006593*t3 - 0.00003169*t4) / 3600);

	// ------------------------------------------------------------------------
	// Compute mean longitude of the mean ascending node of the Moon in radians

	   double om = Math.toRadians((450160.398036 - 6962890.5431*t
	       + 7.4722*t2  + 0.007702*t3 - 0.00005939*t4) / 3600);

	// -----------------------------------------------------------------------
	// Sum 2000B series for nutation in longitude (dPsi) in arc sec * 10000000

	double s = 0;
	s += (-172064161 - 174666*t)*Math.sin(om) + 33386*Math.cos(om);
	s += (-13170906 - 1675*t)*Math.sin(2*(f - d + om)) - 13696*Math.cos(2*(f - d + om));
	s += (-2276413 - 234*t)*Math.sin(2*(f + om)) + 2796*Math.cos(2*(f + om));
	s += (2074554 + 207*t)*Math.sin(2*om) - 698*Math.cos(2*om);
	s += (1475877 - 3633*t)*Math.sin(lp) + 11817*Math.cos(lp);
	s += (-516821 + 1226*t)*Math.sin(lp + 2*(f - d + om)) - 524*Math.cos(lp + 2*(f - d + om));
	s += (711159 + 73*t)*Math.sin(l) - 872*Math.cos(l);
	s += (-387298 - 367*t)*Math.sin(2*f + om) + 380*Math.cos(2*f + om);
	s += (-301461 - 36*t)*Math.sin(l + 2*(f + om)) + 816*Math.cos(l + 2*(f + om));
	s += (215829 - 494*t)*Math.sin(2*(f - d + om) - lp) + 111*Math.cos(2*(f - d + om) - lp);
	s += (128227 + 137*t)*Math.sin(2*(f - d) + om) + 181*Math.cos(2*(f - d) + om);
	s += (123457 + 11*t)*Math.sin(2*(f + om) - l) + 19*Math.cos(2*(f + om) - l);
	s += (156994 + 10*t)*Math.sin(2*d - l) - 168*Math.cos(2*d - l);
	s += (63110 + 63*t)*Math.sin(l + om) + 27*Math.cos(l + om);
	s += (-57976 - 63*t)*Math.sin(om - l) - 189*Math.cos(om - l);
	s += (-59641 - 11*t)*Math.sin(2*(f + d + om) - l) + 149*Math.cos(2*(f + d + om) - l);
	s += (-51613 - 42*t)*Math.sin(l + 2*f + om) + 129*Math.cos(l + 2*f + om);
	s += (45893 + 50*t)*Math.sin(2*(f - l) + om) + 31*Math.cos(2*(f - l) + om);
	s += (63384 + 11*t)*Math.sin(2*d) - 150*Math.cos(2*d);
	s += (-38571 - t)*Math.sin(2*(f + d + om)) + 158*Math.cos(2*(f + d + om));
	s += 32481*Math.sin(2*(f - lp - d + om));
	s -= 47722*Math.sin(2*(d - l)) + 18*Math.cos(2*(d - l));
	s += (-31046 - t)*Math.sin(2*(l + f + om)) + 131*Math.cos(2*(l + f + om));
	s += 28593*Math.sin(l + 2*(f - d + om)) - Math.cos(l + 2*(f - d + om));
	s += (20441 + 21*t)*Math.sin(2*f - l + om) + 10*Math.cos(2*f - l + om);
	s += 29243*Math.sin(2*l) - 74*Math.cos(2*l);
	s += 25887*Math.sin(2*f) - 66*Math.cos(2*f);
	s += (-14053 - 25*t)*Math.sin(lp + om) + 79*Math.cos(lp + om);
	s += (15164 + 10*t)*Math.sin(2*d - l + om) + 11*Math.cos(2*d - l + om);
	s += (-15794 + 72*t)*Math.sin(2*(lp + f - d + om)) - 16*Math.cos(2*(lp + f - d + om));
	s += 21783*Math.sin(2*(d - f)) + 13*Math.cos(2*(d - f));
	s += (-12873 - 10*t)*Math.sin(l - 2*d + om) - 37*Math.cos(l - 2*d + om);
	s += (-12654 + 11*t)*Math.sin(om - lp) + 63*Math.cos(om - lp);
	s -= 10204*Math.sin(2*(f + d) - l + om) - 25*Math.cos(2*(f + d) - l + om);
	s += (16707 - 85*t)*Math.sin(2*lp) - 10*Math.cos(2*lp);
	s -= 7691*Math.sin(l + 2*(f + d + om)) - 44*Math.cos(l + 2*(f + d + om));
	s -= 11024*Math.sin(-2*l + 2*f) + 14*Math.cos(2*(f - l));
	s += (7566 - 21*t)*Math.sin(lp + 2*(f + om)) - 11*Math.cos(lp + 2*(f + om));
	s += (-6637 - 11*t)*Math.sin(2*(f + d) + om) + 25*Math.cos(2*(f + d) + om);
	s += (-7141 + 21*t)*Math.sin(2*(f + om) - lp) + 8*Math.cos(2*(f + om) - lp);
	s += (-6302 - 11*t)*Math.sin(2*d + om) + 2*Math.cos(2*d + om);
	s += (5800 + 10*t)*Math.sin(l + 2*(f - d) + om) + 2*Math.cos(l + 2*(f - d) + om);
	s += 6443*Math.sin(2*(l + f - d + om)) - 7*Math.cos(2*(l + f - d + om));
	s += (-5774 - 11*t)*Math.sin(2*(d - l) + om) - 15*Math.cos(2*(d - l) + om);
	s -= 5350*Math.sin(2*(l + f) + om) - 21*Math.cos(2*(l + f) + om);
	s += (-4752 - 11*t)*Math.sin(2*(f - d) - lp + om) - 3*Math.cos(2*(f - d) - lp + om);
	s += (-4940 - 11*t)*Math.sin(om - 2*d) - 21*Math.cos(om - 2*d);
	s += 7350*Math.sin(2*d - l - lp) - 8*Math.cos(2*d - l - lp);
	s += 4065*Math.sin(2*(l - d) + om) + 6*Math.cos(2*(l - d) + om);
	s += 6579*Math.sin(l + 2*d) - 24*Math.cos(l + 2*d);
	s += 3579*Math.sin(lp + 2*(f - d) + om) + 5*Math.cos(lp + 2*(f - d) + om);
	s += 4725*Math.sin(l - lp) - 6*Math.cos(l - lp);
	s -= 3075*Math.sin(2*(f - l  + om)) + 2*Math.cos(2*(f - l  + om));
	s -= 2904*Math.sin(3*l + 2*(f + om)) - 15*Math.cos(3*l + 2*(f + om));
	s += 4348*Math.sin(2*d - lp) - 10*Math.cos(2*d - lp);
	s -= 2878*Math.sin(l - lp + 2*(f + om)) - 8*Math.cos(l - lp + 2*(f + om));
	s -= 4230*Math.sin(d) - 5*Math.cos(d);
	s -= 2819*Math.sin(2*(f + d + om) - l - lp) - 7*Math.cos(2*(f + d + om) - l - lp);
	s -= 4056*Math.sin(2*f - l) - 5*Math.cos(2*f - l);
	s -= 2647*Math.sin(2*(f + d + om) - lp) - 11*Math.cos(2*(f + d + om) - lp);
	s -= 2294*Math.sin(om - 2*l) + 10*Math.cos(om - 2*l);
	s += 2481*Math.sin(l + lp + 2*(f + om)) - 7*Math.cos(l + lp + 2*(f + om));
	s += 2179*Math.sin(2*l + om) - 2*Math.cos(2*l + om);
	s += 3276*Math.sin(lp - l + d) + Math.cos(lp - l + d);
	s -= 3389*Math.sin(l + lp) - 5*Math.cos(l + lp);
	s += 3339*Math.sin(l + 2*f) - 13*Math.cos(l + 2*f);
	s -= 1987*Math.sin(2*(f - d) - l + om) + 6*Math.cos(2*(f - d) - l + om);
	s -= 1981*Math.sin(l + 2*om);
	s += 4026*Math.sin(d - l) - 353*Math.cos(d - l);
	s += 1660*Math.sin(d + 2*(f + om)) - 5*Math.cos(d + 2*(f + om));
	s -= 1521*Math.sin(2*(f + 2*d + om) - l) - 9*Math.cos(2*(f + 2*d + om) - l);
	s += 1314*Math.sin(lp - l + d + om);
	s -= 1283*Math.sin(2*(f - lp - d) + om);
	s -= 1331*Math.sin(l + 2*(f + d) + om) - 8*Math.cos(l + 2*(f + d) + om);
	s += 1383*Math.sin(2*(f - l + d + om)) - 2*Math.cos(2*(f - l + d + om));
	s += 1405*Math.sin(2*om - l) + 4*Math.cos(2*om - l);
	s += 1290*Math.sin(l + lp + 2*(f - d + om));

	// ---------------------------------------------------
	// Return nutation in ecliptical longitude in degrees.

	double rr = 36000000000.0;
	double dPsiDeg = s / rr;

	return dPsiDeg;

	}
}
