import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

public class HeatMapMain {

	public static void main(String[] args) throws IOException {
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
	}
}
