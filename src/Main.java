import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Scanner;

import javax.imageio.ImageIO;

public class Main {

	// public static final double M_TO_NM = Math.pow(10, -9);
	public static final double J_TO_eV = 6.241509343260179 * Math.pow(10, 18);

	public static final double H_IN_J_S = 6.62607004081 * Math.pow(10, -34);
	public static final double H_IN_eV_S = H_IN_J_S * J_TO_eV;

	public static final double C_IN_M_PER_S = 299792458;

	public static final double K_IN_J_PER_K = 1.3806485279 * Math.pow(10, -23);

	public static final double Q_IN_C = 1.602176620898 * Math.pow(10, -19);

	public static void main(String[] args) throws Exception {
		calculate();
	}

	public static void calculate() {

		double T = 300;
		double vg = 2.7 / H_IN_eV_S;

		double P_L = 1E2;
		double v0 = 2.73 / H_IN_eV_S;
		double a = 1E10;
		double f_L = .95;
		double f = .95;

		double Q_L = Q_L(P_L, v0, vg, a);
		double Q_B = ApproximateBlackBodyIntegral.Q_B(T, vg);

		double V_c = V_c(T);
		double V_op = V_op(V_c, f, Q_L, Q_B);
		double V_g = V_g(vg);

		double z_op = z_op(V_op, V_c);

		double z_m = ApproximateZm.z_m(z_op);
		double m = m(z_m);

		double V_m = V_m(z_m, V_c);

		double ult = ultimateEfficiency(vg, v0, a);
		double nom = nominalEfficiency(T, vg, P_L, v0, a, f_L, f);
		double prac = practicalEfficiency(T, vg, P_L, v0, a, f_L, f);

		System.out.println("Q_L: " + Q_L);
		System.out.println("Q_B: " + Q_B);

		System.out.println("V_c: " + V_c);
		System.out.println("V_g: " + V_g);
		System.out.println("V_op: " + V_op);

		System.out.println("z_op: " + z_op);
		System.out.println("z_m: " + z_m);
		System.out.println("V_m: " + V_m);

		System.out.println("m: " + m);

		System.out.println("ult: " + ult);
		System.out.println("nom: " + nom);
		System.out.println("prac: " + prac);

		// tableOfUltimateEfficiencies(T, vg, P_L, v0, a, f_L, f);
		// tableOfPracticalEfficiencies(T, vg, P_L, v0, a, f_L, f);
		// colorMapUltimateEfficiencies(T, vg, P_L, v0, a, f_L, f);
		// colorMapPracticalEfficiencies(T, vg, P_L, v0, a, f_L, f);

		// graphTemperature(vg, P_L, v0, a, f_L, f);
		// graphLaserPower(T, vg, v0, a, f_L, f);
		// graphFraction(T, vg, P_L, v0, a, f_L);
	}

	public static void graphTemperature(double vg, double P_L, double v0, double a, double f_L, double f) {
		try {
			File outputFile = new File("effVsT.txt");
			outputFile.createNewFile();
			PrintStream outFileStream = new PrintStream(new FileOutputStream(outputFile, false), true);

			String str = "";

			for (int T = 50; T <= 1000; T += 50) {
				double eff = practicalEfficiency(T, vg, P_L, v0, a, f_L, f);
				// System.out.println("T: " + T + " eff: " + eff);
				str += T + "\t" + eff + "\r\n";
			}

			outFileStream.print(str);
			outFileStream.close();
		} catch (Exception e) {

		}
	}

	public static void graphLaserPower(double T, double vg, double v0, double a, double f_L, double f) {
		try {
			File outputFile = new File("effVsP_L.txt");
			outputFile.createNewFile();
			PrintStream outFileStream = new PrintStream(new FileOutputStream(outputFile, false), true);

			String str = "";

			for (double order = 0; order < 10; order += .5) {
				double P_L = Math.pow(10, order);
				double eff = practicalEfficiency(T, vg, P_L, v0, a, f_L, f);
				// System.out.println("P_L: " + P_L + " eff: " + eff);
				str += P_L + "\t" + eff + "\r\n";
			}

			outFileStream.print(str);
			outFileStream.close();
		} catch (Exception e) {

		}
	}

	public static void graphFraction(double T, double vg, double P_L, double v0, double a, double f_L) {
		try {
			File outputFile = new File("effVsf.txt");
			outputFile.createNewFile();
			PrintStream outFileStream = new PrintStream(new FileOutputStream(outputFile, false), true);

			String str = "";

			for (double f = 1E-20; f <= 2; f *= 10) {
				double eff = practicalEfficiency(T, vg, P_L, v0, a, f_L, f);
				// System.out.println("f: " + f + " eff: " + eff);
				str += f + "\t" + eff + "\r\n";
			}

			outFileStream.print(str);
			outFileStream.close();
		} catch (Exception e) {

		}
	}

	public static void tableOfUltimateEfficiencies(double T, double vg, double P_L, double v0, double a, double f_L,
			double f) {
		LaserSemiconductorValues.calculateRealValues("ultimateEfficiencies", new TwoDPercentFunction() {

			@Override
			public double calculate(double v0, double vg) {
				return ultimateEfficiency(vg, v0, a);
			}
		});
	}

	public static void tableOfPracticalEfficiencies(double T, double vg, double P_L, double v0, double a, double f_L,
			double f) {
		LaserSemiconductorValues.calculateRealValues("practicalEfficiencies", new TwoDPercentFunction() {

			@Override
			public double calculate(double v0, double vg) {
				return practicalEfficiency(T, vg, P_L, v0, a, f_L, f);
			}
		});
	}

	public static void colorMapUltimateEfficiencies(double T, double vg, double P_L, double v0, double a, double f_L,
			double f) {

		ColorImageWriter.writeImage(500, 500, "ultEffGrayHeatMap", new TwoDPercentFunction() {

			@Override
			public double calculate(double v0, double vg) {
				double ans = ultimateEfficiency(vg, v0, a);
				return ans;
			}
		});

	}

	public static void colorMapPracticalEfficiencies(double T, double vg, double P_L, double v0, double a, double f_L,
			double f) {

		ColorImageWriter.writeImage(500, 500, "practEffGrayHeatMap", new TwoDPercentFunction() {

			@Override
			public double calculate(double v0, double vg) {
				double ans = practicalEfficiency(T, vg, P_L, v0, a, f_L, f);
				return ans;
			}
		});

	}

	public static double ultimateEfficiency(double vg, double v0, double a) {
		double answer = 1;
		answer *= vg;
		answer /= (Math.PI + 2 * Math.atan(v0 / a));
		answer /= (Math.pow(a, 2) + Math.pow(v0, 2));

		double mainfunc = 0;

		mainfunc += (Math.PI * v0);
		mainfunc += (a * Math.log(1 - 2 * v0 / vg + (Math.pow(a, 2) + Math.pow(v0, 2)) / Math.pow(vg, 2)));
		mainfunc += (2 * v0 * Math.atan((v0 - vg) / a));

		answer *= mainfunc;

		return answer;
	}

	public static double nominalEfficiency(double T, double vg, double P_L, double v0, double a, double f_L, double f) {

		double Q_L = Q_L(P_L, v0, vg, a);
		double Q_B = ApproximateBlackBodyIntegral.Q_B(T, vg);

		double V_c = V_c(T);
		double V_op = V_op(V_c, f, Q_L, Q_B);
		double V_g = V_g(vg);

		double ans = voltageFraction(V_op, V_g) * f_L * ultimateEfficiency(vg, v0, a);
		if (ans < 1)
			return ans;
		else
			return 0;
	}

	public static double practicalEfficiency(double T, double vg, double P_L, double v0, double a, double f_L,
			double f) {
		try {
			if (vg == 0 || v0 == 0)
				return 0;
			double Q_L = Q_L(P_L, v0, vg, a);
			double Q_B = ApproximateBlackBodyIntegral.Q_B(T, vg);

			double V_c = V_c(T);
			double V_op = V_op(V_c, f, Q_L, Q_B);
			double V_g = V_g(vg);

			double z_op = z_op(V_op, V_c);

			double z_m = ApproximateZm.z_m(z_op);

			return m(z_m) * nominalEfficiency(T, vg, P_L, v0, a, f_L, f);
		} catch (Exception e) {
			return 0;
		}
	}

	public static double Q_L(double P_L, double v0, double vg, double a) {
		double ans = 0;
		ans = P_L / (2 * Math.PI * H_IN_J_S);
		ans /= (Math.pow(v0, 2) + Math.pow(a, 2));
		double mainfunc = 0;

		mainfunc += Math.PI * v0;
		mainfunc += a * Math.log((Math.pow(vg - v0, 2) + Math.pow(v0, 2)) / Math.pow(vg, 2));
		mainfunc += +2 * v0 * Math.atan((v0 - vg) / a);

		ans *= mainfunc;
		return ans;
	}

	public static double V_op(double V_c, double f, double Q_L, double Q_B) {
		return V_c * Math.log(f * Q_L / Q_B);
	}

	public static double V_c(double T) {
		double V_c = K_IN_J_PER_K * T / Q_IN_C;
		return V_c;
	}

	public static double V_g(double vg) {
		double V_g = H_IN_J_S * vg / Q_IN_C;
		return V_g;
	}

	public static double z_op(double V_op, double V_c) {
		double z_op = V_op / V_c;
		return z_op;
	}

	public static double V_m(double z_m, double V_c) {
		return z_m * V_c;
	}

	public static double m(double z_m) {
		double ans = Math.pow(z_m, 2);
		ans /= (1 + z_m - Math.exp(-z_m));
		ans /= (z_m + Math.log(1 + z_m));
		return ans;
	}

	public static double voltageFraction(double V_op, double V_g) {
		return V_op / V_g;
	}

	public static String toString(double d) {
		d = round(d, 2);
		return d + "";
	}

	public static double round(double value, int places) {
		if (places < 0)
			throw new IllegalArgumentException();

		BigDecimal bd = new BigDecimal(value);
		bd = bd.setScale(places, RoundingMode.HALF_UP);
		return bd.doubleValue();
	}

}

class ApproximateZm {
	// solves x+ln(1+x)-z_op=0
	// x+ln(1+x)=z_op
	public static double z_m(Double z_op) {
		if (z_op < 0 || z_op.isInfinite() || z_op.isNaN())
			throw new RuntimeException("invalid: " + z_op);
		return solve(z_op, 0, z_op);
	}

	private static double solve(double z_op, double a, double b) {
		double c = (a + b) / 2;
		double c_ans = c + Math.log(c + 1);

		if (percentError(c_ans, z_op) < .0005)
			return c;
		else if (c_ans > z_op)
			return solve(z_op, a, c);
		else
			return solve(z_op, c, b);
	}

	private static double percentError(double est, double val) {
		return Math.abs(est - val) / val;
	}
}

class ApproximateBlackBodyIntegral extends Main {

	public static double Q_B(double T, double vg) {
		double wavelength = C_IN_M_PER_S / vg;
		double ans = integral(T, 0, wavelength, 1000);
		return ans;
	}

	private static double integral(double T, double startWavelength, double endWavelength, double numRect) {
		double rectWidth = endWavelength / numRect;

		double ans = 0;
		for (int i = 0; i < numRect; i++) {
			double xpos1 = i * rectWidth + startWavelength;
			double xpos2 = (i + 1) * rectWidth + startWavelength;
			double midxpos = (xpos1 + xpos2) / 2;

			double val = blackBodyAt(midxpos, T);

			ans += (val * rectWidth);
		}
		return ans;
	}

	private static double blackBodyAt(double wavelength, double T) {
		double exp = Math.exp((H_IN_J_S * C_IN_M_PER_S) / (K_IN_J_PER_K * T * wavelength));
		double ans = 1 / (exp - 1);
		ans /= Math.pow(wavelength, 4);
		ans *= 2 * Math.PI * C_IN_M_PER_S;
		return ans;
	}
}

class ColorImageWriter extends Main {
	public static void writeImage(int w, int h, String name, TwoDPercentFunction function) {
		String path = name + ".png";
		BufferedImage image = new BufferedImage(w, h, BufferedImage.TYPE_BYTE_GRAY);
		double maxeff = 0;
		int maxcol = 0;
		for (int y = 0; y < h; y++) {
			for (int x = 0; x < w; x++) {
				double laserv = ((double) x / w * 8) / H_IN_eV_S;
				double bgv = ((double) y / h * 8) / H_IN_eV_S;
				double efficiency = function.calculate(laserv, bgv);
				if (efficiency > maxeff) {
					maxeff = efficiency;
				}
				int col = (int) (efficiency * 255);
				if (col > maxcol) {
					maxcol = col;
				}
				Color color = new Color(col, col, col);
				image.setRGB(x, y, color.getRGB());
			}
			System.out.println((double) y / h);
		}
		System.out.println(maxeff);
		System.out.println(maxcol);

		File ImageFile = new File(path);
		try {
			ImageIO.write(image, "png", ImageFile);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void imggg() {
		try {
			int width = 500, height = 500;
			int[] pixels1D = new int[width * height];

			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {
					double v0eV = (double) x / width * 8;
					double vgeV = (double) y / height * 8;

					double v0 = v0eV / Main.H_IN_eV_S;
					double vg = vgeV / Main.H_IN_eV_S;

					double efficiency = Main.ultimateEfficiency(vg, v0, Math.pow(10, 9));

					int pixel = (int) (efficiency * 255);

					pixels1D[500 * y + x] = new Color(pixel, pixel, pixel).getRGB();

				}
			}

			BufferedImage pixelImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
			pixelImage.setRGB(0, 0, width, height, pixels1D, 0, width);

			File outFile = new File("test.png");
			outFile.createNewFile();
			ImageIO.write(pixelImage, "png", outFile);

		} catch (Exception e) {
		}
	}

}

abstract class TwoDPercentFunction {

	public abstract double calculate(double a, double b);
}

class LaserSemiconductorValues extends Main {
	public static final File semiconductorFile = new File(
			"C:\\Users\\akars\\Documents\\FIRM\\FIRMProgram\\semiconductors.txt");
	public static final File lasersFile = new File("C:\\Users\\akars\\Documents\\FIRM\\FIRMProgram\\lasers.txt");
	public static final HashMap<String, Double> lasers = new LinkedHashMap<String, Double>();
	public static final HashMap<String, Double> semiconductors = new LinkedHashMap<String, Double>();
	static {

		try {
			Scanner scanSemiconductors = new Scanner(semiconductorFile);
			Scanner scanLasers = new Scanner(lasersFile);

			while (scanLasers.hasNextLine()) {
				String[] strs = scanLasers.nextLine().split(" ");
				try {
					lasers.put(strs[0], Double.parseDouble(strs[1]));
				} catch (Exception e) {
				}
			}

			while (scanSemiconductors.hasNextLine()) {
				String[] strs = scanSemiconductors.nextLine().split(" ");
				try {
					semiconductors.put(strs[0], Double.parseDouble(strs[1]));
				} catch (Exception e) {
				}
			}

			scanSemiconductors.close();
			scanLasers.close();
		} catch (Exception e) {
		}

	}

	public static void calculateRealValues(String fileName, TwoDPercentFunction function) {
		try {
			File outputFile = new File(fileName + ".txt");

			String tabstr = "\t";

			for (String laser : lasers.keySet()) {
				double eV = H_IN_eV_S * C_IN_M_PER_S / lasers.get(laser);

				String eVStr = toString(eV);
				String laserstr = laser + " (" + eVStr + " eV)";

				tabstr += laserstr + "\t";
			}
			tabstr += "\r\n";

			double maxeff = 0;
			String semi = "";
			String laser_ = "";

			for (String semiconductor : semiconductors.keySet()) {

				double bandgap = semiconductors.get(semiconductor);
				double vg = bandgap / H_IN_eV_S;

				String bandgapStr = toString(bandgap);
				tabstr += semiconductor + " (" + bandgapStr + " eV)" + "\t";

				for (String laser : lasers.keySet()) {

					double wavelength = lasers.get(laser);
					double v0 = C_IN_M_PER_S / wavelength;

					double efficiency = function.calculate(v0, vg);

					if (efficiency > maxeff) {
						maxeff = efficiency;
						semi = semiconductor;
						laser_ = laser;
					}
					String effstr = toString(efficiency * 100);

					tabstr += effstr + "%\t";
				}
				tabstr += "\r\n";
			}
			System.out.println(maxeff + " " + laser_ + " " + semi);

			outputFile.createNewFile();
			PrintStream outFileStream = new PrintStream(new FileOutputStream(outputFile, false), true);

			outFileStream.print(tabstr);
			// System.out.println(tabstr);

			outFileStream.close();
		} catch (Exception e) {
		}
	}
}
