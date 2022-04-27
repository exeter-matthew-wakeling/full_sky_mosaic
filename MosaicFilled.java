import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

/**
 * Plot outlines of images that are required to cover a whole sphere mosaic, on a bitmap image.
 * Assume the images are a rectilinear projection, and project the outlines onto either:
 * 1. six rectilinear projections which are the faces of a cube
 * 2. An equirectangular projection
 * Alternatively, produce a list of the image centre coordinates.
 *
 * @author Matthew Wakeling
 */
public class MosaicFilled
{
	public static Random rand = new Random(0);

	public static void main(String[] args) {
		// All angles in degrees.
		double width, height, overlap;
		width = height = overlap = 0.0;
		int type = 0;
		try {
			width = Double.parseDouble(args[0]);
			height = Double.parseDouble(args[1]);
			overlap = Double.parseDouble(args[2]);
			if ("-cross".equals(args[3])) {
				type = 1;
			} else if ("-count".equals(args[3])) {
				type = 2;
			} else if ("-list".equals(args[3])) {
				type = 3;
			} else if (!"-eq".equals(args[3])) {
				throw new RuntimeException("");
			}
		} catch (Exception e) {
			System.err.println("Usage: java MosaicFilled width height overlap type");
			System.err.println("Width and height are the size of the image in degrees - width is the size in the RA direction, and height is the size in the Declination direction.");
			System.err.println("For a standard DSLR camera mounted on a German equatorial mount, width will usually be smaller than height.");
			System.err.println("Overlap is how much the images should overlap in degrees.");
			System.err.println("Type can be:");
			System.err.println("1. -cross - a cross-rectilinear projection pnm image will be written to the stdout.");
			System.err.println("2. -eq - an equirectangular projection pnm image will be written to the stdout.");
			System.err.println("3. -count - reports the number of images required.");
			System.err.println("4. -list - lists the coordinates of the image centres.");
			System.exit(1);
		}
		// Output image is a 1600x1200 bitmap for cross-rectilinear, and 1800x900 for equirectangular. 4 values necessary - RGB and image count.
		int[] image = new int[type == 0 ? 1800 * 900 * 4 : 1600 * 1200 * 4];
		// Arrangement: Make one image pointed directly at each pole. Then have a set of ring layers going between them.
		// The number of images in each ring will depend on the declination of that ring.
		double ringSpace = 180.0 - Math.min(width, height);
		// Number of image layers required in the declination direction, from pole to pole:
		int decCount = (int) Math.ceil(ringSpace / (height - overlap));
		double decOverlap = (height * decCount - ringSpace) / (decCount + 1);
		double[] sagAboves = new double[decCount];
		double[] sagBelows = new double[decCount];
		boolean needMoreRings = true;
		while (needMoreRings) {
			// The exact spacing of the rings depends on how much sag there is in each ring. The sag in each ring depends on the exact spacing of the rings.
			// So, we need to do a few iterations until these converge on a solution.
			for (int i = 0; i < 7; i++) {
				double totalSag = 0.0;
				//System.err.println("Trying " + decCount + " rings plus the polar images. Dec overlap: " + decOverlap);
				double startDec = -ringSpace / 2.0 - decOverlap;
				int[] raCounts = new int[decCount];
				for (int dec = 0; dec < decCount; dec++) {
					double decA = startDec + height / 2.0 - sagBelows[dec];
					startDec += height - sagBelows[dec] - sagAboves[dec] - decOverlap;
					// How wide in RA is each image at its narrowest point, after projection:
					double minWidth = minimalWidth(width - overlap, height, decA);
					// Number of images in the layer:
					int raCount = (int) Math.ceil(360.0 / minWidth);
					raCounts[dec] = raCount;
					double decSag = decA - height / 2.0;
					// Sag is the declination difference between decSag and the position of the top/bottom edge of the image where it overlaps with its neighbour.
					// It is decSag minus the declination of the equator rotated around the x axis by decSag at ra 180/raCount.
					// Original position of equator: x = sin(ra), y = 0, z = cos(ra)
					// Rotated position of equator: x = sin(ra), y = cos(ra) * sin(decSag), z = cos(ra) * cos(decSag)
					// Find value of atan(y, sqrt(x*x+z*z)) when x/z = tan(180/raCount) = sin(ra) / (cos(ra) * cos(decSag)). tan(ra) = tan(180/raCount) * cos(decSag)
					double raSagRadians = Math.atan(Math.tan(Math.PI / raCount) * Math.cos(Math.PI * decSag / 180.0));
					//double sag = Math.abs(decSag) * (1.0 - Math.cos(Math.PI / raCount));
					double sagBelow = decSag < 0.0 ? 180.0 * Math.asin(Math.cos(raSagRadians) * Math.sin(Math.PI * decSag / 180.0)) / Math.PI - decSag : 0.0;
					decSag = decA + height / 2.0;
					raSagRadians = Math.atan(Math.tan(Math.PI / raCount) * Math.cos(Math.PI * decSag / 180.0));
					double sagAbove = decSag > 0.0 ? decSag - 180.0 * Math.asin(Math.cos(raSagRadians) * Math.sin(Math.PI * decSag / 180.0)) / Math.PI : 0.0;
					sagBelows[dec] = sagBelow;
					sagAboves[dec] = sagAbove;
					//System.err.println("Ring " + dec + " dec " + decA + " minWidth " + minWidth + " count " + raCount + " spacing " + (360.0 / raCount) + " sagBelow " + sagBelow + " sagAbove " + sagAbove);
					totalSag += sagBelow + sagAbove;
				}
				if (raCounts[0] == 4) {
					// The polar image has four edges too, so we can remove the sag.
					totalSag -= sagBelows[0];
					sagBelows[0] = 0.0;
				}
				if (raCounts[decCount - 1] == 4) {
					totalSag -= sagAboves[decCount - 1];
					sagAboves[decCount - 1] = 0.0;
				}
				for (int dec = 0; dec < decCount - 1; dec++) {
					if (raCounts[dec] == raCounts[dec + 1]) {
						// Two rings have the same number of images, we don't need the sag between them
						totalSag -= sagAboves[dec];
						totalSag -= sagBelows[dec + 1];
						sagAboves[dec] = 0.0;
						sagBelows[dec + 1] = 0.0;
					}
				}
				decOverlap = (height * decCount - ringSpace - totalSag) / (decCount + 1);
				//System.err.println("Total sag: " + totalSag + ", actual dec overlap: " + decOverlap);
			}
			needMoreRings = false;
			if (decOverlap < overlap) {
				// If the total amount of sag is too much to allow enough overlap in the dec direction, then we need to add another ring.
				decCount++;
				sagBelows = new double[decCount];
				sagAboves = new double[decCount];
				decOverlap = (height * decCount - ringSpace) / (decCount + 1);
				needMoreRings = true;
			}
		}
		System.err.println("Using " + decCount + " rings plus the polar images. Overlap between rings is " + decOverlap);
		int totalImages = 2;
		double startDec = -ringSpace / 2.0 - decOverlap;
		for (int dec = 0; dec < decCount; dec++) {
			double decA = startDec + height / 2.0 - sagBelows[dec];
			startDec += height - sagBelows[dec] - sagAboves[dec] - decOverlap;
			double minWidth = minimalWidth(width - overlap, height, decA);
			int raCount = (int) Math.ceil(360.0 / minWidth);
			totalImages += raCount;
			System.err.println("Ring " + dec + " declination " + decA + " count " + raCount + " overlap " + (minWidth + overlap - 360.0 / raCount));
		}
		System.err.println("Total " + totalImages + " images required");
		printOutline(0.0, -90.0, width, height, image, type);
		startDec = -ringSpace / 2.0 - decOverlap;
		for (int dec = 0; dec < decCount; dec++) {
			double decA = startDec + height / 2.0 - sagBelows[dec];
			startDec += height - sagBelows[dec] - sagAboves[dec] - decOverlap;
			// Number of images in the layer:
			double minWidth = minimalWidth(width - overlap, height, decA);
			int raCount = (int) Math.ceil(360.0 / minWidth);
			for (int ra = 0; ra < raCount; ra++) {
				double raA = (ra * 360.0) / raCount;
				printOutline(raA, decA, width, height, image, type);
			}
		}
		printOutline(0.0, 90.0, width, height, image, type);

		if (type == 0) {
			System.out.println("P6 1800 900 255");
		} else if (type == 1) {
			System.out.println("P6 1600 1200 255");
		}
		if ((type == 0) || (type == 1)) {
			int xsize = type == 0 ? 1800 : 1600;
			for (int y = 0; y < (type == 0 ? 900 : 1200); y++) {
				for (int x = 0; x < xsize; x++) {
					double c = image[(x + y * xsize) * 4 + 3];
					c = c < 2 ? 1.5 : c;
					int r = c == 0 ? 0 : (int) (image[(x + y * xsize) * 4] / c);
					int g = c == 0 ? 0 : (int) (image[(x + y * xsize) * 4 + 1] / c);
					int b = c == 0 ? 0 : (int) (image[(x + y * xsize) * 4 + 2] / c);
					System.out.write(r);
					System.out.write(g);
					System.out.write(b);
				}
			}
			System.out.flush();
		}
		if (type == 2) {
			System.out.println(totalImages);
		}
	}

	public static double minimalWidth(double width, double height, double dec) {
		if (dec == 0.0) {
			return width;
		}
		dec = Math.abs(dec);
		if (dec <= height / 2.0) {
			// The image overlaps the equator, so the narrowest point is the corner of the image.
			return 360.0 * Math.atan(Math.tan(width * Math.PI / 360.0) / (Math.cos(dec * Math.PI / 180.0) + Math.tan(height * Math.PI / 360.0) * Math.sin(dec * Math.PI / 180.0))) / Math.PI;
		}
		// Find the width of the image in ra at its narrowest point
		// The side edge of the image is a line with x = tan(width / 2), y = -tan(h), z = 1, where h is the height position in the image, down from the centre.
		// Rotated by dec, this is x = tan(width / 2), y = sin(dec) - tan(h) * cos(dec), z = cos(dec) + tan(h) * sin(dec)
		// Declination of the point is atan(y/sqrt(x*x+z*z)) = atan((sin(dec) - tan(h) * cos(dec)) / sqrt(tan(width / 2) * tan(width / 2) + (cos(dec) + tan(h) * sin(dec)) * (cos(dec) + tan(h) * sin(dec))))
		// We need to search for a value of h where this declination is dec - height / 2, and a first approximation for this is height / 2.
		// Iteratively refine the value.
		double sinDec = Math.sin(dec * Math.PI / 180.0);
		double cosDec = Math.cos(dec * Math.PI / 180.0);
		double tanWid = Math.tan(width * Math.PI / 360.0);
		double target = dec - height / 2.0;
		double h = height / 2.0;
		double lastH = 0.0;
		while (Math.abs(h - lastH) > 0.0000001) {
			double tanH = Math.tan(h * Math.PI / 180.0);
			double tDec = 180 * Math.atan((sinDec - tanH * cosDec) / Math.sqrt(tanWid * tanWid + (cosDec + tanH * sinDec) * (cosDec + tanH * sinDec))) / Math.PI;
			lastH = h;
			h += tDec - target;
		}
		// Now that we have our h, we can calculate the width as 2 * atan(x / z) = 2 * atan(tan(width/2) / (cos(dec) + tan(h) * sin*dec))
		double tanH = Math.tan(h * Math.PI / 180.0);
		return 360.0 * Math.atan(tanWid / (cosDec + tanH * sinDec)) / Math.PI;
	}

	public static void printOutline(double ra, double dec, double width, double height, int[] image, int type) {
		if (type == 3) {
			System.out.println(ra + "\t" + dec + "\t" + ((int) (ra / 15.0)) + "h " + (((int) (ra * 4.0)) % 60) + "' " + (((int) (ra * 240.0)) % 60) + "\"\t" + ((int) dec) + "\u00b0 " + (((int) Math.abs(dec * 60.0)) % 60) + "' " + (((int) Math.abs(dec * 3600.0)) % 60) + "\"");
			return;
		}
		if (type == 2) {
			// Just counting them
			return;
		}
		double tanWidth = Math.tan(width * Math.PI / 360.0);
		double tanHeight = Math.tan(height * Math.PI / 360.0);
		double cosDec = Math.cos(dec * Math.PI / 180.0);
		double sinDec = Math.sin(dec * Math.PI / 180.0);
		double firstRa = 0.0;
		double firstDec = 0.0;
		double lastRa = 0.0;
		double lastDec = 0.0;
		// Draw a rectangle, projected, as a set of lines, stored in this collection:
		ArrayList<Line> lines = new ArrayList<Line>();
		for (double i = -tanWidth; i <= tanWidth; i += 0.01) {
			double x = i;
			double y = tanHeight;
			double z = 1.0;
			// Rotate around x axis by dec.
			double newY = cosDec * y + sinDec * z;
			double newZ = cosDec * z - sinDec * y;
			double newDec = 180.0 * Math.atan2(newY, Math.sqrt(x * x + newZ * newZ)) / Math.PI;
			double newRa = ra + 180.0 * Math.atan2(x, newZ) / Math.PI;
			while (newRa < -180.0) {
				newRa += 360.0;
			}
			while (newRa >= 180.0) {
				newRa -= 360.0;
			}
			if (i == -tanWidth) {
				firstRa = newRa;
				firstDec = newDec;
			} else {
				lines.add(new Line(lastRa, lastDec, newRa, newDec));
			}
			lastRa = newRa;
			lastDec = newDec;
		}
		for (double i = tanHeight; i >= -tanHeight; i -= 0.01) {
			double x = tanWidth;
			double y = i;
			double z = 1.0;
			// Rotate around x axis by dec.
			double newY = cosDec * y + sinDec * z;
			double newZ = cosDec * z - sinDec * y;
			double newDec = 180.0 * Math.atan2(newY, Math.sqrt(x * x + newZ * newZ)) / Math.PI;
			double newRa = ra + 180.0 * Math.atan2(x, newZ) / Math.PI;
			while (newRa < -180.0) {
				newRa += 360.0;
			}
			while (newRa >= 180.0) {
				newRa -= 360.0;
			}
			lines.add(new Line(lastRa, lastDec, newRa, newDec));
			lastRa = newRa;
			lastDec = newDec;
		}
		for (double i = tanWidth; i >= -tanWidth; i -= 0.01) {
			double x = i;
			double y = -tanHeight;
			double z = 1.0;
			// Rotate around x axis by dec.
			double newY = cosDec * y + sinDec * z;
			double newZ = cosDec * z - sinDec * y;
			double newDec = 180.0 * Math.atan2(newY, Math.sqrt(x * x + newZ * newZ)) / Math.PI;
			double newRa = ra + 180.0 * Math.atan2(x, newZ) / Math.PI;
			while (newRa < -180.0) {
				newRa += 360.0;
			}
			while (newRa >= 180.0) {
				newRa -= 360.0;
			}
			lines.add(new Line(lastRa, lastDec, newRa, newDec));
			lastRa = newRa;
			lastDec = newDec;
		}
		for (double i = -tanHeight; i <= tanHeight; i += 0.01) {
			double x = -tanWidth;
			double y = i;
			double z = 1.0;
			// Rotate around x axis by dec.
			double newY = cosDec * y + sinDec * z;
			double newZ = cosDec * z - sinDec * y;
			double newDec = 180.0 * Math.atan2(newY, Math.sqrt(x * x + newZ * newZ)) / Math.PI;
			double newRa = ra + 180.0 * Math.atan2(x, newZ) / Math.PI;
			while (newRa < -180.0) {
				newRa += 360.0;
			}
			while (newRa >= 180.0) {
				newRa -= 360.0;
			}
			lines.add(new Line(lastRa, lastDec, newRa, newDec));
			lastRa = newRa;
			lastDec = newDec;
		}
		lines.add(new Line(lastRa, lastDec, firstRa, firstDec));

		// Now we have a list of lines, but some of them wrap around the -180 to 180 RA discontinuity. Whether this matters or not depends on the projection.

		// Now we have a list of lines. Scanline implementation of filled polygon rendering.
		int r = 50 + rand.nextInt(206);
		int g = 50 + rand.nextInt(206);
		int b = 50 + rand.nextInt(206);
		//System.err.println("R: " + r + ", G: " + g + ", B: " + b);
		if (type == 0) {
			// Equirectangular projection, 1800x900 resolution.
			// RA wraparound matters, so we need to deal with it.
			ArrayList<Line> newLines = new ArrayList<Line>();
			for (Line line : lines) {
				if ((line.x1 < -90.0) && (line.x2 > 90.0)) {
					// Line points leftwards through the meridian - duplicate it
					newLines.add(new Line(line.x1, line.y1, line.x2 - 360.0, line.y2));
					newLines.add(new Line(line.x1 + 360.0, line.y1, line.x2, line.y2));
				} else if ((line.x1 > 90.0) && (line.x2 < -90)) {
					// Line points rightwards through the meridian - duplicate it
					newLines.add(new Line(line.x1 - 360.0, line.y1, line.x2, line.y2));
					newLines.add(new Line(line.x1, line.y1, line.x2 + 360, line.y2));
				} else {
					newLines.add(line);
				}
			}
			lines = newLines;
			// There is some complication caused by the top and bottom images, so it is best to do scanning up-down.
			for (int x = 0; x < 1800; x++) {
				// Work out which lines intersect with this x value.
				double xRa = x * 0.2 - 180.0;
				ArrayList<Double> intersections = new ArrayList<Double>();
				for (Line line : lines) {
					if (((line.x1 <= xRa) && (line.x2 > xRa)) || ((line.x1 > xRa) && (line.x2 <= xRa))) {
						double inter = line.y1 + (xRa - line.x1) * (line.y2 - line.y1) / (line.x2 - line.x1);
						intersections.add(inter);
					}
				}
				if (!intersections.isEmpty()) {
					Collections.sort(intersections);
					//System.err.println("Intersections: " + intersections);
					if (dec < 0.0) {
						// Scan from positive dec downwards
						for (int i = intersections.size() - 1; i >= 0; i -= 2) {
							double startDec = intersections.get(i);
							double endDec = i == 0 ? -90.0 : intersections.get(i - 1);
							int startY = (int) (450.0 - 5.0 * startDec + 0.5);
							int endY = (int) (450.0 - 5.0 * endDec - 0.5);
							//System.err.println("Starty: " + startY + ", endY: " + endY + ", x: " + x);
							for (int y = startY; y <= endY; y++) {
								image[(x + y * 1800) * 4] += r;
								image[(x + y * 1800) * 4 + 1] += g;
								image[(x + y * 1800) * 4 + 2] += b;
								image[(x + y * 1800) * 4 + 3]++;
							}
						}
					} else {
						// Scan from negative dec upwards
						for (int i = 0; i < intersections.size(); i += 2) {
							double startDec = intersections.get(i);
							double endDec = i == intersections.size() - 1 ? 90.0 : intersections.get(i + 1);
							int startY = (int) (450.0 - 5.0 * endDec + 0.5);
							int endY = (int) (450.0 - 5.0 * startDec - 0.5);
							//System.err.println("Starty: " + startY + ", endY: " + endY + ", x: " + x);
							for (int y = startY; y <= endY; y++) {
								image[(x + y * 1800) * 4] += r;
								image[(x + y * 1800) * 4 + 1] += g;
								image[(x + y * 1800) * 4 + 2] += b;
								image[(x + y * 1800) * 4 + 3]++;
							}
						}
					}
				}
			}
		} else if (type == 1) {
			// Cross-rectilinear projection.
			// Need to make 6 separate rectilinear projections, each 400x400, on a 1600x1200 image.
			crossRectProject(image, 0, 400, 400, lines, r, g, b);
			crossRectProject(image, 1, 800, 400, lines, r, g, b);
			crossRectProject(image, 2, 1200, 400, lines, r, g, b);
			crossRectProject(image, 3, 0, 400, lines, r, g, b);
			crossRectProject(image, 4, 400, 0, lines, r, g, b);
			crossRectProject(image, 5, 400, 800, lines, r, g, b);
		}
	}

	public static void crossRectProject(int[] image, int rotation, int xoff, int yoff, ArrayList<Line> lines, int r, int g, int b) {
		ArrayList<Line> projected = new ArrayList<Line>();
		boolean isVisible = false;
		for (Line line : lines) {
			double[] proj1 = crossRectProjectPoint(rotation, line.x1, line.y1);
			double[] proj2 = crossRectProjectPoint(rotation, line.x2, line.y2);
			if ((proj1[2] > 0.5) || (proj2[2] > 0.5)) {
				isVisible = true;
			}
			if ((proj1[2] > 0.0) && (proj2[2] > 0.0)) {
				double x1 = 200.0 + 200.0 * proj1[0] / proj1[2];
				double y1 = 200.0 - 200.0 * proj1[1] / proj1[2];
				double x2 = 200.0 + 200.0 * proj2[0] / proj2[2];
				double y2 = 200.0 - 200.0 * proj2[1] / proj2[2];
				projected.add(new Line(x1, y1, x2, y2));
			}
		}
		if (isVisible) {
			for (int y = 0; y < 400; y++) {
				// Work out which lines intersect with this y value.
				ArrayList<Double> intersections = new ArrayList<Double>();
				Set<Double> positiveInter = new HashSet<Double>();
				for (Line line : projected) {
					if (((line.y1 <= y) && (line.y2 > y)) || ((line.y1 > y) && (line.y2 <= y))) {
						double inter = line.x1 + (y - line.y1) * (line.x2 - line.x1) / (line.y2 - line.y1);
						intersections.add(inter);
						if (line.y1 > line.y2) {
							positiveInter.add(inter);
						}
					}
				}
				if (!intersections.isEmpty()) {
					Collections.sort(intersections);
					int startX = 0;
					boolean inPolygon = false;
					if (!positiveInter.contains(intersections.get(0))) {
						inPolygon = true;
					}
					for (int i = 0; i < intersections.size(); i++) {
						if (positiveInter.contains(intersections.get(i))) {
							inPolygon = true;
							startX = Math.max(0, Math.min(400, (int) (intersections.get(i) + 0.5)));
						} else {
							inPolygon = false;
							int endX = Math.min(400, (int) (intersections.get(i) + 0.5));
							for (int x = startX; x < endX; x++) {
								image[(x + xoff + (y + yoff) * 1600) * 4] += r;
								image[(x + xoff + (y + yoff) * 1600) * 4 + 1] += g;
								image[(x + xoff + (y + yoff) * 1600) * 4 + 2] += b;
								image[(x + xoff + (y + yoff) * 1600) * 4 + 3]++;
							}
						}
					}
					if (inPolygon) {
						for (int x = startX; x < 400; x++) {
							image[(x + xoff + (y + yoff) * 1600) * 4] += r;
							image[(x + xoff + (y + yoff) * 1600) * 4 + 1] += g;
							image[(x + xoff + (y + yoff) * 1600) * 4 + 2] += b;
							image[(x + xoff + (y + yoff) * 1600) * 4 + 3]++;
						}
					}
				}
			}
		}
	}

	public static double[] crossRectProjectPoint(int rotation, double ra, double dec) {
		double x = 0.0;
		double y = 0.0;
		double z = 0.0;
		switch (rotation) {
			case 0:
				x = Math.cos(Math.PI * dec / 180.0) * Math.sin(Math.PI * ra / 180.0);
				y = Math.sin(Math.PI * dec / 180.0);
				z = Math.cos(Math.PI * dec / 180.0) * Math.cos(Math.PI * ra / 180.0);
				break;
			case 1:
				x = Math.cos(Math.PI * dec / 180.0) * Math.sin(Math.PI * (ra - 90.0) / 180.0);
				y = Math.sin(Math.PI * dec / 180.0);
				z = Math.cos(Math.PI * dec / 180.0) * Math.cos(Math.PI * (ra - 90.0) / 180.0);
				break;
			case 2:
				x = Math.cos(Math.PI * dec / 180.0) * Math.sin(Math.PI * (ra - 180.0) / 180.0);
				y = Math.sin(Math.PI * dec / 180.0);
				z = Math.cos(Math.PI * dec / 180.0) * Math.cos(Math.PI * (ra - 180.0) / 180.0);
				break;
			case 3:
				x = Math.cos(Math.PI * dec / 180.0) * Math.sin(Math.PI * (ra + 90.0) / 180.0);
				y = Math.sin(Math.PI * dec / 180.0);
				z = Math.cos(Math.PI * dec / 180.0) * Math.cos(Math.PI * (ra + 90.0) / 180.0);
				break;
			case 4:
				x = Math.cos(Math.PI * dec / 180.0) * Math.sin(Math.PI * ra / 180.0);
				y = -Math.cos(Math.PI * dec / 180.0) * Math.cos(Math.PI * ra / 180.0);
				z = Math.sin(Math.PI * dec / 180.0);
				break;
			case 5:
				x = Math.cos(Math.PI * dec / 180.0) * Math.sin(Math.PI * ra / 180.0);
				y = Math.cos(Math.PI * dec / 180.0) * Math.cos(Math.PI * ra / 180.0);
				z = -Math.sin(Math.PI * dec / 180.0);
				break;
		}
		return new double[] {x, y, z};
	}

	public static class Line
	{
		public double x1, y1, x2, y2;

		public Line(double x1, double y1, double x2, double y2) {
			this.x1 = x1;
			this.y1 = y1;
			this.x2 = x2;
			this.y2 = y2;
		}
	}
}
