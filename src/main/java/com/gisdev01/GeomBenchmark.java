package com.gisdev01;

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.geotools.data.DataStore;
import org.geotools.data.DataStoreFinder;
import org.geotools.data.FeatureSource;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.CoordinateList;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LinearRing;
import com.vividsolutions.jts.geom.MultiPolygon;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.Polygon;
import com.vividsolutions.jts.geom.prep.PreparedPolygon;
import com.vividsolutions.jts.io.ParseException;
import com.vividsolutions.jts.io.WKTReader;
import com.vividsolutions.jts.io.WKTWriter;
import com.vividsolutions.jts.simplify.DouglasPeuckerSimplifier;

// Note: We use fully qualified class names from com.esri.core.geometry due to several name collisions with the above
// JTS names (eg. Polygon, Point, etc.)


public class GeomBenchmark {

	// test shapefile can be downloaded here:
	// https://www.census.gov/cgi-bin/geo/shapefiles/index.php (as of Oct. 2018)
	private static String SHAPEFILE_COUNTIES_FILEPATH = "test/resources/tl_2018_us_county.shp";

	private List<Polygon> jtsPolygons;
	private List<Envelope> jtsEnvelopes;

	private List<Envelope> jtsClipBoxes;
	private List<Polygon> ellipses;
	private List<Integer> ids;

	private List<PreparedPolygon> ppolygons;
	
	private List<com.esri.core.geometry.Polygon> esriPolygons;
	private List<com.esri.core.geometry.Envelope> esriEnvelopes;

	private List<com.esri.core.geometry.Envelope2D> esriEnvelopes2DClipBoxes; // for clips
	private List<com.esri.core.geometry.Polygon> esriPolygonEllipses; // for intersections
	private List<Integer> esriIds;

	private List<com.esri.core.geometry.Polygon> esriPreparedPolys;

    private com.esri.core.geometry.OperatorFactoryLocal factory = com.esri.core.geometry.OperatorFactoryLocal.getInstance();
    private com.esri.core.geometry.OperatorImportFromWkt operatorImport = (com.esri.core.geometry.OperatorImportFromWkt) factory.getOperator(com.esri.core.geometry.Operator.Type.ImportFromWkt);
    private com.esri.core.geometry.OperatorExportToWkt operatorExport = (com.esri.core.geometry.OperatorExportToWkt) factory.getOperator(com.esri.core.geometry.Operator.Type.ExportToWkt);
    private com.esri.core.geometry.OperatorIntersection operatorIntersection = (com.esri.core.geometry.OperatorIntersection) factory.getOperator(com.esri.core.geometry.Operator.Type.Intersection);
    private com.esri.core.geometry.OperatorConvexHull operatorConvexHull = (com.esri.core.geometry.OperatorConvexHull) factory.getOperator(com.esri.core.geometry.Operator.Type.ConvexHull);
    private com.esri.core.geometry.OperatorClip operatorClip = (com.esri.core.geometry.OperatorClip) factory.getOperator(com.esri.core.geometry.Operator.Type.Clip);
    private com.esri.core.geometry.OperatorSimplify operatorSimplify = (com.esri.core.geometry.OperatorSimplify) factory.getOperator(com.esri.core.geometry.Operator.Type.Simplify);
    private com.esri.core.geometry.OperatorWithin operatorWithin = (com.esri.core.geometry.OperatorWithin) factory.getOperator(com.esri.core.geometry.Operator.Type.Within);
    private com.esri.core.geometry.OperatorContains operatorContains = (com.esri.core.geometry.OperatorContains) factory.getOperator(com.esri.core.geometry.Operator.Type.Contains);
    private com.esri.core.geometry.OperatorGeneralize operatorGeneralize = (com.esri.core.geometry.OperatorGeneralize) factory.getOperator(com.esri.core.geometry.Operator.Type.Generalize);
    private com.esri.core.geometry.SpatialReference sr = com.esri.core.geometry.SpatialReference.create(4269);
	private WKTWriter wkt = new WKTWriter();
	private WKTReader wktr = new WKTReader();
	private GeometryFactory geometryFactory;


    public static void main(String[] args) throws Exception {
        GeomBenchmark geomBenchmark = new GeomBenchmark();
        File shapefileFile = new File(SHAPEFILE_COUNTIES_FILEPATH);
        geomBenchmark.prepare(shapefileFile);
        geomBenchmark.runchecks();
    }

	private GeomBenchmark() {
		super();

		// jts objects
		jtsPolygons = new ArrayList<Polygon>();
		jtsEnvelopes = new ArrayList<Envelope>();
		jtsClipBoxes = new ArrayList<Envelope>();
		ellipses = new ArrayList<Polygon>();
		geometryFactory = new GeometryFactory();
		ids = new ArrayList<Integer>();

		// esri objects
		esriPolygons = new ArrayList<com.esri.core.geometry.Polygon>();
		esriEnvelopes = new ArrayList<com.esri.core.geometry.Envelope>();
		esriEnvelopes2DClipBoxes = new ArrayList<com.esri.core.geometry.Envelope2D>();
		esriPolygonEllipses = new ArrayList<com.esri.core.geometry.Polygon>();
		esriIds = new ArrayList<Integer>();

	}

	
	private com.esri.core.geometry.Polygon JTStoESRIPoly(Polygon jtsGeom){
		String wktP = wkt.write(jtsGeom);
		return (com.esri.core.geometry.Polygon) operatorImport.execute(0, com.esri.core.geometry.Geometry.Type.Polygon , wktP, null);
	}
	
	public void prepare(File shapefile) {
		jtsPolygons.clear();
		jtsEnvelopes.clear();
		jtsClipBoxes.clear();
		ellipses.clear();
		ids.clear();
		
		esriPolygons.clear();
		esriEnvelopes.clear();
		esriEnvelopes2DClipBoxes.clear();
		esriPolygonEllipses.clear();
		esriIds.clear();

		geometryFactory = new GeometryFactory();

		try {
			// load all shapes from shapefile into jtsPolygons
			readShapefile(shapefile, jtsPolygons, ids);

		} catch (Exception e) {
			System.out.println(e.getMessage());
			System.exit(1);
		}

		ppolygons = new ArrayList<PreparedPolygon>(jtsPolygons.size());
		

		for (Polygon polygon : jtsPolygons) {
			PreparedPolygon prepoly = new PreparedPolygon(polygon);
			ppolygons.add(prepoly);
			if (prepoly.contains(polygon.getCentroid())) {
				// needed to cache it
			}
		}

		// Create envelopes
		for (Polygon polygon : jtsPolygons) {
			jtsEnvelopes.add(polygon.getEnvelopeInternal());
			com.esri.core.geometry.Polygon p2 = JTStoESRIPoly(polygon);
			esriPolygons.add(p2);
			com.esri.core.geometry.Envelope e2 = new com.esri.core.geometry.Envelope();
			p2.queryEnvelope(e2);
			esriEnvelopes.add(e2);
		}
		
		
		// Create the star-ellipses for intersections later on
		if (Compare.MEASURE_OVERLAY || Compare.MEASURE_CLIP) {
			int k = 0;
			for (Envelope box : jtsEnvelopes) {
				k++;

				double cx = box.centre().x;
				double cy = box.centre().y;

				double dx = box.getWidth();
				double dy = box.getHeight();

				if (Compare.MEASURE_OVERLAY) {
					double a1 = Compare.OVERLAY_ELLIPSE_FACTOR1 * 0.5 * dx;
					double b1 = Compare.OVERLAY_ELLIPSE_FACTOR1 * 0.5 * dy;
					double a2 = Compare.OVERLAY_ELLIPSE_FACTOR2 * 0.5 * dx;
					double b2 = Compare.OVERLAY_ELLIPSE_FACTOR2 * 0.5 * dy;

					// We will use a coordinate list to build the linearring
					CoordinateList clist = new CoordinateList();
					// Compare.OVERLAY_ELLIPSE_COUNT);
					double angle = 0.0; // 45.0 * ggl::math::d2r; //0.0;
					for (int i = 0; i < Compare.OVERLAY_ELLIPSE_COUNT - 1; i++, angle += Compare.delta) {
						if (i % 2 == 0) {
							clist.add(new Coordinate(cx + a1 * Math.sin(angle),
									cy + b1 * Math.cos(angle)));
						} else {
							clist.add(new Coordinate(cx + a2 * Math.sin(angle),
									cy + b2 * Math.cos(angle)));
						}
					}

					clist.add(clist.get(0));
					LinearRing lr = geometryFactory.createLinearRing(clist.toCoordinateArray());
					Polygon ellipse = geometryFactory.createPolygon(lr, null);
					ellipses.add(ellipse);
					com.esri.core.geometry.Polygon e2 = JTStoESRIPoly(ellipse);
					esriPolygonEllipses.add(e2);
				}

				if (Compare.MEASURE_CLIP) {
					// note : weird use of sin/cos for a constant angle
					// effectively this create a box . shrinkBy(
					// Compare.CLIP_FACTOR * 0.5*sqrt2)?
					double a = Compare.CLIP_FACTOR * 0.5 * dx;
					double b = Compare.CLIP_FACTOR * 0.5 * dy;

					double angle1 = Math.toRadians(225.0);
					double angle2 = Math.toRadians(45.0);

					double x0 = (cx + a * Math.sin(angle1));
					double y0 = (cy + b * Math.cos(angle1));

					double x1 = (cx + a * Math.sin(angle2));
					double y1 = (cy + b * Math.cos(angle2));

					Envelope clipbox = new Envelope(x0, x1, y0, y1);
					jtsClipBoxes.add(clipbox);
					com.esri.core.geometry.Envelope2D env = new com.esri.core.geometry.Envelope2D(x0, y0, x1, y1);
					esriEnvelopes2DClipBoxes.add(env);
				}
			}
		}
	}

	
	public void runchecks() {
		if (Compare.MEASURE_AREA) {
			System.out.println("");
			//JTS
			double area = 0;
			long t0 = System.nanoTime();
			for (int i = 0; i < Compare.AREA_COUNT; i++) {
				for (Polygon polygon : jtsPolygons) {
					area += polygon.getArea();
				}
			}
			long t1 = System.nanoTime();
			Compare.report_area("JTS", t1 - t0, jtsPolygons.size(), area);
			
			//ESRI
			area = 0;
			t0 = System.nanoTime();
			for (int i = 0; i < Compare.AREA_COUNT; i++) {
				for (com.esri.core.geometry.Polygon polygon : esriPolygons) {
					area += polygon.calculateArea2D();
				}
			}
			t1 = System.nanoTime();
			Compare.report_area("ESRI", t1 - t0, jtsPolygons.size(), area);
		}

		if (Compare.MEASURE_CENTROID) {
			System.out.println("");
			//JTS
			double sum_x = 0, sum_y = 0;
			long t0 = System.nanoTime();
			for (int i = 0; i < Compare.CENTROID_COUNT; i++) {
				for (Polygon polygon : jtsPolygons) {
					Point centroid = polygon.getCentroid();
					sum_x += centroid.getX();
					sum_y += centroid.getY();
				}
			}
			long t1 = System.nanoTime();
			Compare.report_centroid("JTS", t1 - t0, jtsPolygons.size(), sum_x, sum_y);
			
			//ESRI
			Compare.report_centroid("ESRI - UNSUPPORTED", -1L, -1, -1D, -1D);
			
		}

		if (Compare.MEASURE_CONVEX_HULL) {
			System.out.println("");
			//JTS
			double area = 0.0;
			long t0 = System.nanoTime();
			for (Polygon polygon : jtsPolygons) {
				Geometry hull = polygon.convexHull();
				if (Compare.HULL_AREA) {
					area += Math.abs(hull.getArea());
				}
			}
			long t1 = System.nanoTime();
			Compare.report_hull("JTS", t1 - t0, jtsPolygons.size(), area);
			
			//ESRI
			area = 0.0;
			t0 = System.nanoTime();
			for (com.esri.core.geometry.Polygon polygon : esriPolygons) {
				com.esri.core.geometry.Polygon hull = (com.esri.core.geometry.Polygon) operatorConvexHull.execute(polygon,null);
				if (Compare.HULL_AREA) {
					area += Math.abs(hull.calculateArea2D());
				}
			}
			t1 = System.nanoTime();
			Compare.report_hull("ESRI", t1 - t0, jtsPolygons.size(), area);
		}

		if (Compare.MEASURE_OVERLAY) {
			System.out.println("");
			//JTS
			double area1 = 0.0, area2 = 0.0;
			long t0 = System.nanoTime();
			for (int i = 0; i < Compare.OVERLAY_COUNT; i++) {
				int k = 0;
				Iterator<Polygon> eit = ellipses.iterator();
				for (Iterator<Polygon> pit = jtsPolygons.iterator(); pit.hasNext()	&& eit.hasNext(); k++) {
					Polygon poly = pit.next();
					Polygon ellipse = eit.next();
					if (Compare.OVERLAY_AREA) {
						area1 += poly.getArea();
					}
					Geometry v = ellipse.intersection(poly);
					if (Compare.OVERLAY_AREA) {
						area2 += v.getArea();
					}
				}
			}
			long t1 = System.nanoTime();
			Compare.report_overlay("JTS", t1 - t0, jtsPolygons.size(), area1, area2);
			
			
			//ESRI
			t0 = System.nanoTime();
			area1 = 0;
			area2 = 0;
			for (int i = 0; i < Compare.OVERLAY_COUNT; i++) {
				int k = 0;
				Iterator<com.esri.core.geometry.Polygon> eit = esriPolygonEllipses.iterator();
				for (Iterator<com.esri.core.geometry.Polygon> pit = esriPolygons.iterator(); pit.hasNext()&& eit.hasNext(); k++) {
					com.esri.core.geometry.Polygon poly = pit.next();
					com.esri.core.geometry.Polygon ellipse = eit.next();
					if (Compare.OVERLAY_AREA) {
						area1 += poly.calculateArea2D();
					}
					com.esri.core.geometry.GeometryCursor cursor1 = new com.esri.core.geometry.SimpleGeometryCursor(poly);
					com.esri.core.geometry.GeometryCursor cursor2 = new com.esri.core.geometry.SimpleGeometryCursor(ellipse);
					com.esri.core.geometry.GeometryCursor outputGeoms = operatorIntersection.execute(cursor1, cursor2, sr, null);
					com.esri.core.geometry.Geometry intersect = outputGeoms.next();
					if (Compare.OVERLAY_AREA) {
						area2 += intersect.calculateArea2D();
					}
				}
			}
			t1 = System.nanoTime();
			Compare.report_overlay("ESRI", t1 - t0, jtsPolygons.size(), area1, area2);
			
		}

		if (Compare.MEASURE_CLIP) {
			System.out.println("");
			//JTS
			boolean first = true;
			double area1 = 0.0, area2 = 0.0;
			long t0 = System.nanoTime();
			for (int i = 0; i < Compare.CLIP_COUNT; i++) {
				Iterator<Envelope> bit = jtsClipBoxes.iterator();
				Iterator<Polygon> pit = jtsPolygons.iterator();
				for (int k = 0; pit.hasNext() && bit.hasNext(); k++) {
					Polygon poly = pit.next();
					Envelope clipenv = bit.next();
					Geometry clipgeom = geometryFactory.toGeometry(clipenv);
					if (Compare.CLIP_AREA) {
						area1 += poly.getArea();
					}
					Geometry v = clipgeom.intersection(poly);
					if (Compare.CLIP_AREA) {
						area2 += v.getArea();
					}
				}
			}
			long t1 = System.nanoTime();
			Compare.report_clip("JTS", t1 - t0, jtsPolygons.size(), area1, area2);
			
			//ESRI
			first = true;
			area1 = 0.0; 
			area2 = 0.0;
			t0 = System.nanoTime();
			for (int i = 0; i < Compare.CLIP_COUNT; i++) {
				Iterator<com.esri.core.geometry.Envelope2D> bit = esriEnvelopes2DClipBoxes.iterator();
				Iterator<com.esri.core.geometry.Polygon> pit = esriPolygons.iterator();
				for (int k = 0; pit.hasNext() && bit.hasNext(); k++) {
					com.esri.core.geometry.Polygon poly = pit.next();
					com.esri.core.geometry.Envelope2D clipenv = bit.next();

					
					if (Compare.CLIP_AREA) {
						area1 += poly.calculateArea2D();
					}
					com.esri.core.geometry.Geometry g = operatorClip.execute(poly, clipenv, sr, null);
					//com.esri.core.geometry.Geometry v = operatorIntersection.execute(poly,g, sr, null);
					if (Compare.CLIP_AREA) {
						area2 += g.calculateArea2D();
					}
				}
			}
			t1 = System.nanoTime();
			Compare.report_clip("ESRI", t1 - t0, jtsPolygons.size(), area1, area2);
		}

		if (Compare.MEASURE_SIMPLIFY) {
			System.out.println("");
			//JTS
			int count1 = 0, count2 = 0;
			double length1 = 0.0, length2 = 0.0;
			long t0 = System.nanoTime();
			for (Polygon polygon : jtsPolygons) {
				Geometry simplegeom = DouglasPeuckerSimplifier.simplify(
						polygon, Compare.SIMPLIFY_DISTANCE);
				count1 += polygon.getNumPoints();
				count2 += simplegeom.getNumPoints();
				if (Compare.SIMPLIFY_LENGTH) {
					length1 += polygon.getLength();
					length2 += simplegeom.getLength();
				}

			}
			long t1 = System.nanoTime();
			Compare.report_simplify("JTS", t1 - t0, jtsPolygons.size(), length1, length2,	count1, count2);
			
			//ESRI
			count1 = 0;
			count2 = 0;
			length1 = 0.0;
			length2 = 0.0;
			t0 = System.nanoTime();
			for (com.esri.core.geometry.Polygon polygon : esriPolygons) {
				
				//com.esri.core.geometry.Geometry simplegeom = operatorSimplify.execute(polygon, sr, false,null);
				com.esri.core.geometry.Geometry simplegeom = operatorGeneralize.execute(polygon, Compare.SIMPLIFY_DISTANCE, false, null);
				count1 += polygon.getPointCount();
				count2 += ((com.esri.core.geometry.Polygon)simplegeom).getPointCount();
				if (Compare.SIMPLIFY_LENGTH) {
					length1 += polygon.calculateLength2D();
					length2 += simplegeom.calculateLength2D();
				}

			}
			t1 = System.nanoTime();
			Compare.report_simplify("ESRI", t1 - t0, jtsPolygons.size(), length1, length2,	count1, count2);

		}

		
		if (Compare.MEASURE_WITHIN) {
			System.out.println("");
			//JTS
			int count = 0;
			long t0 = System.nanoTime();
			for (int e = 0; e < jtsEnvelopes.size(); e++) {
				Envelope b = jtsEnvelopes.get(e);
				Coordinate c = b.centre();
				Point p = geometryFactory.createPoint(c);
				Iterator<Envelope> bit = jtsEnvelopes.iterator();
				Iterator<Polygon> pit = jtsPolygons.iterator();
				for (int k = 0; pit.hasNext() && bit.hasNext(); k++) {
					Polygon poly = pit.next();
					Envelope box = bit.next();
					if (box.contains(c) && p.within(poly)) {
						count++;
					}
				}
			}
			long t1 = System.nanoTime();
			Compare.report_within("JTS", t1 - t0, jtsPolygons.size(), count, -1);
			
			
			//ESRI
			count = 0;
			t0 = System.nanoTime();
			for (int e = 0; e < esriEnvelopes.size(); e++) {
				com.esri.core.geometry.Envelope b = esriEnvelopes.get(e);
				
				com.esri.core.geometry.Point p = b.getCenter();
				Iterator<com.esri.core.geometry.Envelope> bit = esriEnvelopes.iterator();
				Iterator<com.esri.core.geometry.Polygon> pit = esriPolygons.iterator();
				for (int k = 0; pit.hasNext() && bit.hasNext(); k++) {
					com.esri.core.geometry.Polygon poly = pit.next();
					com.esri.core.geometry.Envelope box = bit.next();
					if (box.contains(p) && operatorWithin.execute(p, poly, sr,null)) {
						count++;
					}
				}
			}
			t1 = System.nanoTime();
			Compare.report_within("ESRI", t1 - t0, jtsPolygons.size(), count, -1);
		}


		if (Compare.MEASURE_CONTAINS) {
			System.out.println("");
			//JTS
			int count = 0;
			List<Point> points = new ArrayList<Point>(jtsEnvelopes.size());
			for (int e = 0; e < jtsEnvelopes.size(); e++) {
				Envelope b = jtsEnvelopes.get(e);
				Coordinate c = b.centre();
				Point p = geometryFactory.createPoint(c);
				points.add(p);
			}
			long t0 = System.nanoTime();
			Iterator<Point> pointIt = points.iterator();
			Iterator<Polygon> pit = jtsPolygons.iterator();
			for (int k = 0; pit.hasNext() && pointIt.hasNext(); k++) {
				Polygon poly = pit.next();
				Point p = pointIt.next();
				if (poly.contains(p)) {
					count++;
				}
			}
			long t1 = System.nanoTime();
			Compare.report_contains("JTS", t1 - t0, jtsPolygons.size(), count, -1);
			
			
			
			
			//ESRI
			count = 0;
			List<com.esri.core.geometry.Point> Epoints = new ArrayList<com.esri.core.geometry.Point>(jtsEnvelopes.size());
			for (int e = 0; e < esriEnvelopes.size(); e++) {
				com.esri.core.geometry.Envelope b = esriEnvelopes.get(e);
				com.esri.core.geometry.Point p = b.getCenter();
				Epoints.add(p);
			}
			
			t0 = System.nanoTime();
			Iterator<com.esri.core.geometry.Point> EpointIt = Epoints.iterator();
			Iterator<com.esri.core.geometry.Polygon> Epit = esriPolygons.iterator();
			for (int k = 0; Epit.hasNext() && EpointIt.hasNext(); k++) {
				com.esri.core.geometry.Polygon poly = Epit.next();
				com.esri.core.geometry.Point p = EpointIt.next();
				if (operatorContains.execute(poly, p, sr, null)) {
					count++;
				}
			}
			t1 = System.nanoTime();
			Compare.report_contains("ESRI", t1 - t0, jtsPolygons.size(), count, -1);

		}
		
		//robustness check
		//http://tsusiatsoftware.net/jts/jts-faq/jts-faq.html#D
		//Intersection of LINESTRING(0 0, 5 3), LINESTRING(0 0, 1.2 0.72) should be (0 0, 1.2, 0.72)
		//JTS
		String wkt1 = "LINESTRING(0 0, 5 3 )";
		String wkt2 = "LINESTRING(0 0, 1.2 0.72)";
		
		Geometry gi1 = null;
		Geometry gi2 = null;
		try {
			gi1 = wktr.read(wkt1);
		    gi2 = wktr.read(wkt2);
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		System.out.println("");
		System.out.println("TEST of the intersction of LINESTRING(0 0, 5 3) & LINESTRING(0 0, 1.2 0.72)");
		System.out.println("Actual value is LINESTRING(0 0, 1.2 0.72) - but issues occur due to finite precision - checks robustness");
		System.out.println("JTS Intersection: " + wkt.writeFormatted(gi1.intersection(gi2)));
		
		//ESRI
		com.esri.core.geometry.Geometry Egi1 = (com.esri.core.geometry.Geometry) operatorImport.execute(0, com.esri.core.geometry.Geometry.Type.Unknown , wkt1, null);
		com.esri.core.geometry.Geometry Egi2 = (com.esri.core.geometry.Geometry) operatorImport.execute(0, com.esri.core.geometry.Geometry.Type.Unknown , wkt2, null);
		com.esri.core.geometry.GeometryCursor c1 = new com.esri.core.geometry.SimpleGeometryCursor(Egi1);
		com.esri.core.geometry.GeometryCursor c2 = new com.esri.core.geometry.SimpleGeometryCursor(Egi2);
		com.esri.core.geometry.Geometry inter12 = operatorIntersection.execute(c1, c2, sr, null).next();
		System.out.println("ESRI Intersection: " + operatorExport.execute(0, inter12, null));
		
		
	}

	private static void readShapefile(File file, List<Polygon> polygons,
                                      List<Integer> ids) {
		try {
			/*
			 * Attempt to find a GeoTools DataStore that can handle the shapefile
			 */
			Map<String, Serializable> connectParameters = new HashMap<String, Serializable>();

			connectParameters.put("url", file.toURI().toURL());
			connectParameters.put("create spatial index", false);

			DataStore dataStore = DataStoreFinder
					.getDataStore(connectParameters);
			if (dataStore == null) {
				System.out.println("No DataStore found to handle" + file.getPath());
				System.exit(1);
			}

			/*
			 * We are now connected to the shapefile. Get the type name of the
			 * features within it
			 */
			String[] typeNames = dataStore.getTypeNames();
			String typeName = typeNames[0];

			System.out.println("Reading content " + typeName);

			/*
			 * Iterate through the features, collecting some spatial data (line
			 * or boundary length) on each one
			 */
			FeatureSource<SimpleFeatureType, SimpleFeature> featureSource;
			FeatureCollection<SimpleFeatureType, SimpleFeature> collection;
			FeatureIterator<SimpleFeature> iterator;

			featureSource = dataStore.getFeatureSource(typeName);
			collection = featureSource.getFeatures();
			iterator = collection.features();

			int i = 0;
			double totalArea = 0.0;
			try {
				while (iterator.hasNext()) {
					SimpleFeature feature = iterator.next();

					i++;
					// System.out.println(i+"  :"+feature.getID());
					/*
					 * The spatial portion of the feature is represented by a
					 * Geometry object
					 */
					Geometry geometry = (Geometry) feature.getDefaultGeometry();
					// TODO validate polygon?

					// Process only jtsPolygons, and from them only single-jtsPolygons
					// without holes
					Polygon polygon = null;
					if (geometry instanceof Polygon) {
						polygon = (Polygon) geometry;
					} else if (geometry instanceof MultiPolygon) {
						MultiPolygon mp = (MultiPolygon) geometry;
						if (mp.getNumGeometries() == 1) {
							polygon = (Polygon) mp.getGeometryN(0);
						} else {
							/*System.out
									.println(i
											+ "skipped:  not a single polygon multipolygon");
							*/
						}
					} else {
						System.out.println(i
								+ " skipped: not a (multi)polygon:"
								+ geometry.getGeometryType());
					}

					if (polygon != null) {
						if (polygon.getNumInteriorRing() == 0) {
							totalArea += polygon.getArea();
							polygons.add(polygon);
							ids.add(i);
						} else {
							/*System.out.println(i
									+ "  not a single ring polygon:"
									+ geometry.getGeometryType());
							*/
						}
					}

				}
			} finally {
				/*
				 * You MUST explicitly close the feature iterator otherwise
				 * terrible things will happen !!!
				 */
				if (iterator != null) {
					iterator.close();
				}
			}

			System.out.println("Total Area: " + totalArea);

		} catch (Exception ex) {
			ex.printStackTrace();
			System.exit(1);
		}
	}
}
