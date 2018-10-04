Originally based on this previous work: https://github.com/chrisbennight/jts-esri-benchmark

## To test locally:
Download one of the U.S. Census Tiger Counties Shapefiles (multiple files make up a shapefile (.shp)). 
https://www.census.gov/cgi-bin/geo/shapefiles/index.php
Select a Layer Type: Counties
Place the unzipped contents of the downloaded .zip file into the ./test/resources/ directory.

## To Run:
<pre>
mvn clean compile
mvn exec:exec
</pre>


## Output

<pre>
[INFO] Scanning for projects...
[INFO]                                                                         
[INFO] ------------------------------------------------------------------------
[INFO] Building geom-performance-check 0.1-SNAPSHOT
[INFO] ------------------------------------------------------------------------
[INFO] 
[INFO] --- exec-maven-plugin:1.2.1:exec (default-cli) @ geom-performance-check ---
Reading content c_02ap13
Total Area730.0445917511591

JTS AREA: 2803 total: 730.0445917511679 other: 0 time: 0.073104899 s
ESRI AREA: 2803 total: 730.0445917511674 other: 0 time: 0.08710406 s

JTS CENTROID: 2803 total: -91.60642950541535 other: 38.08915874601412 time: 0.277397862 s
ESRI - UNSUPPORTED CENTROID: -1 total: 0.1 other: 0.1 time: -1.0E-9 s

JTS HULL: 2803 total: 822.4889226894915 other: 0 time: 0.578264693 s
ESRI HULL: 2803 total: 822.4889226894913 other: 0 time: 1.394710517 s

JTS INTERSECTION: 101 total: 158.5542576517064 other: 730.0445917511591 time: 6.606177194 s
ESRI INTERSECTION: 101 total: 158.55425765345527 other: 730.0445917511589 time: 25.704174385 s

JTS CLIP: 2803 total: 370.1493045493158 other: 730.0445917511591 time: 0.945118263 s
ESRI CLIP: 2803 total: 370.14930454931573 other: 730.0445917511589 time: 0.514534584 s

JTS SIMPLIFY: 2803 total: 6197.906981556354 other: 6213.17050681722 time: 0.752725677 s
ESRI SIMPLIFY: 2803 total: 6198.08245123076 other: 6213.17050681722 time: 0.357620154 s

JTS WITHIN: 2803 total: 2796 other: -1 time: 0.843608158 s
ESRI WITHIN: 2803 total: 2796 other: -1 time: 0.454011268 s

JTS CONTAINS: 2803 total: 2761 other: -1 time: 0.252709395 s
ESRI CONTAINS: 2803 total: 2761 other: -1 time: 0.062269236 s

TEST of the intersction of LINESTRING(0 0, 5 3) & LINESTRING(0 0, 1.2 0.72)
Actual value is LINESTRING(0 0, 1.2 0.72) - but issues occur due to finite precision - checks robustness
JTS Intersection: LINESTRING (0 0, 1.2 0.72)
ESRI Intersection: MULTILINESTRING ((0 0, 1.2000000000000002 0.72))
[INFO] ------------------------------------------------------------------------
[INFO] BUILD SUCCESS
[INFO] ------------------------------------------------------------------------
[INFO] Total time: 52.720s
[INFO] Finished at: Sun Jun 09 17:47:57 EDT 2013
[INFO] Final Memory: 6M/183M
[INFO] ------------------------------------------------------------------------
</pre>
