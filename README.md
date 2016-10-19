# netstat
Network K-funcition calculation

Steps of work:
1.	compile the source code to .exe
2.	place all .shp, .txt, .exe files in one directory with net_stat_cmd.exe
3.	split roads by Sanet(or ArcGIS), if not sure in the quality of road file
4.	in cmd change directory to the one with all the files
5.	convert shapes to txt by shpLIB in cmd:
"shpdump.exe shpname.shp > txtname.txt"
6.	run the programm in cmd:
"net_stat_cmd.exe" and follow requests
7.	build the plots in your software of choice (e.g. Excel)
