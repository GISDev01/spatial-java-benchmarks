package com.gisdev01;

public class Compare {

	public final static int AREA_COUNT = 10;

	public final static double AREA_DIVISOR = 1.0;
	public final static double LENGTH_DIVISOR = 1.0;

	public final static double SIMPLIFY_DISTANCE = 0.0008975;
	// results in simplified length of 99%, vertex reduction of ca. 10%

	public final static int CENTROID_COUNT = 10;

	public final static double CLIP_DISTANCE = 0.1;

	public final static int OVERLAY_COUNT = 1;
	public final static int CLIP_COUNT = 1;
	public final static double INTEGER_FACTOR = 10000.0;

	// #define OVERLAY_UNION

	public final static int OVERLAY_ELLIPSE_COUNT = 101;// 101;
	public final static double delta = Math.PI * 2.0/(OVERLAY_ELLIPSE_COUNT - 1);
	public final static double OVERLAY_ELLIPSE_FACTOR1 = 1.1; // 1.1
	public final static double OVERLAY_ELLIPSE_FACTOR2 = 0.2; // 0.2
	public final static double CLIP_FACTOR = 0.9;

	public final static boolean MEASURE_AREA = true;
	public final static boolean MEASURE_CENTROID = true;
	public final static boolean MEASURE_CLIP = true;
	public final static boolean MEASURE_CONVEX_HULL = true;
	public final static boolean MEASURE_OVERLAY = true; //effe uitezet ivm 
	public final static boolean MEASURE_SIMPLIFY = true;
	public final static boolean MEASURE_TOUCH = true;
//	 currently only for GEOS, not further worked out

	public final static boolean MEASURE_CONTAINS = true; // jts ONLY so far
	public final static boolean MEASURE_WITHIN = true;

	// Variables to check the results
	public final static boolean CLIP_AREA = true;
	public final static boolean HULL_AREA = true;
	public final static boolean OVERLAY_AREA = true;
	public final static boolean SIMPLIFY_LENGTH = true;

	public static void report_area(String type, long l, int n, double area) {
		area /= AREA_COUNT;
		area /= AREA_DIVISOR;
		report_common(l, n, type + " AREA:", area, 0);
	}

	public static void report_centroid(String type, long l, int n, double sumX, double sumY) {
		sumX /= n * CENTROID_COUNT;
		sumY /= n * CENTROID_COUNT;
		report_common(l, n, type + " CENTROID:", sumX, sumY);
	}

	public static void report_hull(String type, long l, int n, double area) {
		report_common(l, n, type + " HULL:", area / AREA_DIVISOR, 0);
	}

	public static void report_clip(String type, long l, int n, double area1, double area2) {
		area1 /= CLIP_COUNT;
		area2 /= CLIP_COUNT;
		area1 /= AREA_DIVISOR;
		area2 /= AREA_DIVISOR;
		report_common(l, n, type + " CLIP:", area2, area1);
	}

	public static void report_overlay(String type, long l, int size, double area1, double area2) {
		area1 /= OVERLAY_COUNT;
		area2 /= OVERLAY_COUNT;
		area1 /= AREA_DIVISOR;
		area2 /= AREA_DIVISOR;
		report_common(l, OVERLAY_ELLIPSE_COUNT,
		// #if defined(OVERLAY_UNION)
				// "UNION"
				// #elif defined(OVERLAY_SYMDIFF)
				// "SYM_DIFFERENCE"
				// #else
				type + " INTERSECTION:"
				// #endif
				, area2, area1);
	}

	public static void report_simplify(String type, long l, int n, double length1,	double length2, int count1, int count2) {
		length1 /= LENGTH_DIVISOR;
		length2 /= LENGTH_DIVISOR;
		report_common(l, n, type + " SIMPLIFY:", length2, length1);
	}

	public static void report_within(String type, long l, int size, int count, int count_boundary) {
		report_common(l, size, type + " WITHIN:", count, count_boundary);
	}

	public static void report_contains(String type, long l, int size, int count,	int count_boundary) {
		report_common(l, size, type + " CONTAINS:", count, count_boundary);
	}

	public static void report_touch(String type, long l, int n, int count, int count_box) {
		report_common(l, n, type + " TOUCH:", count, count_box);
	}

	public static void report_common(double s, int n, String algorithm, Object value1, Object value2) {
		System.out.println(algorithm + " " + n + " total: " + value1+ " other: " + value2 + " time: " + s / 1.0e9 + " s");
		System.gc();
		try {
			Thread.currentThread().sleep(100);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
