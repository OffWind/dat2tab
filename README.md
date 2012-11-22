dat2tab
=======

Convert a wind time series data file (.DAT) to Observed Wind Climate file (.TAB)

Usage:

> dat2tab <filename.dat>

Output: filename.dat.tab

Input file format:

< Header line 1: with met station name, coordinates, etc. Content is indifferent.>
< Header line 2: Period of time series. Content is indifferent.>
< Header line 3: Number of records in data series. Content is indifferent.>
< Header line 4: Column headings. Content is indifferent.>
< Data lines 5 to (nrecords+4). Wind data order is mandatory and must be comma-separated:

col 1-station name, 
col 2-mean wind direction in 10min window, 
col 3-mean wind velocity in 10min window, 
col 4-max. speed in 10min window
col 5-min. speed in 10min window
col 6-std. dev. of wind speed in 10min window
col 7-time stamp

Example DAT file:

MAST 1 LON=-8.54512, LAT=42.312876, H=80m agl (WGS84)
January to December 2006
records  52560
Mast, dir, Vave, Vmax, Vmin, stdev, YYYYMMDDhhmm,
44,310,11.5,18.3,5.3,2.9,200601010000,
44,315,14,20.2,6.5,2.9,200601010010,
44,310,13,18.7,6.1,2.7,200601010020,
44,313,11.3,18.7,1.9,2.9,200601010030,
...