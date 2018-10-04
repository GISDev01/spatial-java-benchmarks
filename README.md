Testing performance on various workloads across various Java spatial libraries. 

Compiled/tested against JDK 8 Update 152

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


## Sample console output with Census 2018 Counties data
Placeholder.

Originally based on this previous work: https://github.com/chrisbennight/jts-esri-benchmark
