#include <stdio.h> /* standard input-output */
#include <stdlib.h> /* for memory allocation */
#include <math.h> /* for sin, arctg... */
#include <ctype.h> /* for isdigit */


#define max(A,B) ((A) > (B) ? (A) : (B)) /* macroses for extent estimation */
#define min(A,B) ((A) < (B) ? (A) : (B))


/* declare external variables */

int a, cr_num = 0, rd_num = 0, nd_num, cr_null;
int *w, *k, *p, *p1; /* arrays of result values */
double *b, *kb, *pk, *det, *dt;
char *pointfile , *linefile; /* input text files of points and lines */

double tot_len, st = 50, d = 50, h, max_dist = 5000;
double **cr_dist, **cr_euc, **cr_det; /* network and planar distance matrices of points */

FILE *fp, *fp1;

struct POINT {
	double x, y;
	};

struct CRASH {
	struct POINT pt; 
	double r_dist, pf_dist, pl_dist;
	struct ROAD *adj_road;
	struct CRASH *more;	/* to place several crashes on the same segment, then on the whole link */
	struct NODE *pf, *pl;
	} *crashes;
	
struct ROAD {
	struct POINT ext1, ext2, pf, pl;
	struct ROAD *link;
	struct CRASH *crashes;
	} *roads;

struct NODE {
	struct POINT pt;
	struct LINK *link;
	double tot_dist;
	} *nodes;

struct LINK {
	struct NODE *node;
	double length;
	struct CRASH *crash;
	struct LINK *more;
	};

int get_crashes();
int get_roads();
int crashes_to_roads();
int adj_list();
int dist_matrix();
int result();
int crash_segment(int, int, struct ROAD *);
int get_coord(FILE*, double*, double*);
struct ROAD *add_segment(FILE*, struct ROAD *);
double dist(struct POINT, struct POINT);
int go(struct NODE *, int);
double pt_dist(struct POINT, struct ROAD);

int main (int argc, char *argv[])
{
	int i, j;

	struct ROAD *segm;
	struct CRASH *r;
	struct LINK *s;

	printf("Print the point txt file name:\n");

	scanf("%s", pointfile);

	// pointfile = argv[1];

	printf("Print the line txt file name, m:\n");

	scanf("%s", linefile);

	// linefile = argv[2];

	printf("Print the starting search distance, m:\n");

	scanf("%lf", &st);

	// st = strtod (argv[3], NULL);

	printf("Print the increment of the search distance, m:\n");

	scanf("%lf", &d);

	// d = strtod (argv[4], NULL);

	printf("Print the maximum search distance, m:\n");

	scanf("%lf", &max_dist);

	// max_dist = strtod (argv[5], NULL);

//	FILE *outp;

//	outp = fopen("rez.txt", "w");

	printf ("Reading point file..\n");

	get_crashes();

	printf ("Reading line file..\n");

	get_roads();

	printf ("Assigning points to lines..\n");

	crashes_to_roads();

	printf ("Building graph..\n");

	adj_list();


	
	printf ("Calculating distances..\n");

	dist_matrix();

	printf ("Writing results..\n");

	result();
	
	return 0;
}


/* read the .txt file of the crashes to the array of CRASH structures */

int get_crashes()
{
	int i=0;

	/* read, how many objects */

	fp = fopen(pointfile, "r");
		
	while ( (a=fgetc(fp)) != '#' )
		;
	while ((a=fgetc(fp)) != '\n' )
		if isdigit(a)
			cr_num = cr_num * 10 + (a - '0'); 

	/* allocate memory */

	crashes = malloc(sizeof(struct CRASH) * cr_num);

	/* read the file */
	
	while ( (a=fgetc(fp)) != EOF ) 
	{
		/* skip all white and unnecessary symbols */
		if (a == '\n')
		{	for (a = fgetc(fp); a == ' ' || a == '	' ; a=fgetc(fp))
				;
			/* check if the string starts with '(' */
			if (a == '(')
			{	get_coord(fp, &crashes[i].pt.x, &crashes[i].pt.y);
				crashes[i].r_dist = 70.0; /* in meters - buffer to search */
				crashes[i].adj_road = NULL; 
				crashes[i].more = NULL;
				i++;
			}					
		}
	}

	fclose (fp);

	return 0;

}


/* read the .txt file of the roads to the array of ROAD structures */

int get_roads()
{
	int i=0, n=1;

	/* read, how many objects, then allocate memory */

	fp = fopen(linefile, "r");
		
	while ( (a=fgetc(fp)) != '#' )
		;
	while ((a=fgetc(fp)) != '\n' )
		if isdigit(a)
			rd_num = rd_num * 10 + (a - '0'); 

	roads = malloc(sizeof(struct ROAD) * rd_num);

	/* read the file */
      
	while ( (a=fgetc(fp)) != EOF ) 
	{
		/* skip all white and unnecessary symbols until the end of line */
		if (a == '\n')
		{		/* skip all white and unnecessary symbols in the begining of line, assign first non-white symbol to a */
				for (a = fgetc(fp); a == ' ' || a == '	' ; a=fgetc(fp))
		
			/* check what the string starts with */
			if (a == 'B')
				get_coord(fp, &roads[i].ext1.x , &roads[i].ext1.y );				
			if (a == 't')
				get_coord(fp, &roads[i].ext2.x , &roads[i].ext2.y );
			
			if ((a == '(') && (n == 1))
			{	get_coord(fp, &roads[i].pf.x , &roads[i].pf.y );
				n++;
			}
			if ((a == '(') && (n == 2))
			{	get_coord(fp, &roads[i].pl.x , &roads[i].pl.y );
				n++;
			} 
			roads[i].link = NULL;
			if ((a == '(') && (n == 3))
			{	roads[i].link = add_segment(fp, roads+i);
			}

			if ((a == '\n' || a == EOF) && n != 1)
			{	
				roads[i].crashes = NULL;
				i++;
				n=1;
			}
		}
	}

	fclose(fp);

	return 0;

}


/* assign the crashes to the roads */

int crashes_to_roads()
{
	int i, j;
	struct ROAD *segm;

	/* for every road */
	for (i=0; i < rd_num; i++)
				
		for (j=0; j < cr_num; j++)
		
			/* check every crash if it is within its extent+70m */

			if (crashes[j].pt.x <= max (roads[i].ext1.x, roads[i].ext2.x) + 70.0 && crashes
			[j].pt.x >= min (roads[i].ext1.x, roads[i].ext2.x) - 70.0 && crashes[j].pt.y <= 
			max (roads[i].ext1.y, roads[i].ext2.y) + 70.0 && crashes[j].pt.y >= min (roads
			[i].ext1.y, roads[i].ext2.y) - 70.0)

			{

				/* check all the links */
				for (segm = (roads + i); segm != NULL; segm = segm->link)
					crash_segment(i, j, segm);

			}
/* check how many crashes are with NULL adj_roads */
	cr_null=0;
	for (j=0; j < cr_num; j++)
		if (crashes[j].adj_road == NULL)
			cr_null++;
	printf ("-- %d points are not assigned to lines..\n", cr_null);

	return 0;

}


/* Build a linked list of nodes with the pointers to the crashes and distances */

int adj_list()
{

	int i, j, k=1, n = 0;
	double len;

	struct POINT f, l;
	struct ROAD *segm;
	struct CRASH *cr, *r;
	struct LINK **g;

	/* from the number of 'roads' objects learn preliminary, maximum size of 'nodes' array */
	
	nodes = malloc(sizeof(struct NODE) * 2 * rd_num);

	/* for every road */
	for (i=0; i < rd_num; i++)
	{
		/* create the two pt - f, l; length len, *crash cr for the link */
		len = 0.0;
		f = roads[i].pf;
		r = cr = NULL;
		/* go until the last vertex of the road */
		for (segm = roads + i; segm != NULL; segm = segm->link)
		{			
			/* check for crashes - put them with distance to the first vertex */
			if (segm->crashes != NULL)
			{	/* add the crashes of the segment to the crashes of the road */
				if (cr == NULL)
					r = cr = segm->crashes;
				else	r->more = segm->crashes;
				
				while (r->more != NULL)
				{	r->pf_dist = len + dist(segm->pf, r->pt); 
					r = r->more;
				}
				r->pf_dist = len + dist(segm->pf, r->pt); 
			}
			/* calculate the distance */
			len = len + dist(segm->pf, segm->pl);
			l = segm->pl;
		}

		tot_len += len;

		for (j=0; j != n; j++)
		{	if (nodes[j].pt.x == f.x && nodes[j].pt.y == f.y)
			{	for (k=j+1; k != n; k++)
					if (nodes[k].pt.x == l.x && nodes[k].pt.y == l.y)
						break;
				if (k == n) 
				{	nodes[k].pt.x = l.x;
					nodes[k].pt.y = l.y;
					nodes[k].link = NULL;
					n++;
				}
				break;
			}

			if (nodes[j].pt.x == l.x && nodes[j].pt.y == l.y)
			{	for (k=j+1; k != n; k++)
					if (nodes[k].pt.x == f.x && nodes[k].pt.y == f.y)
						break;
				if (k == n) 
				{	nodes[k].pt.x = f.x;
					nodes[k].pt.y = f.y;
					nodes[k].link = NULL;
					n++;
				}
				break;
			}
		}

		if (j == n)
		{	nodes[j].pt.x = f.x;
			nodes[j].pt.y = f.y;
			nodes[j].link = NULL;
			n++;
			nodes[k=j+1].pt.x = l.x;
			nodes[k].pt.y = l.y;
			nodes[k].link = NULL;
			n++;
 		}


		/* finish creating the crashes */
		if (nodes[j].pt.x == f.x && nodes[j].pt.y == f.y)
			for (r = cr; r != NULL; r = r->more)
			{	r->pf = nodes+j;
				r->pl = nodes+k;
				r->pl_dist = len - r->pf_dist;
			}
		else
			for (r = cr; r != NULL; r = r->more)
			{	r->pf = nodes+k;	
				r->pl = nodes+j;
				r->pl_dist = len - r->pf_dist;
			}

		/* add the link */
		for (g = &nodes[j].link; *g != NULL; g = &(*g)->more)
			;
		*g = malloc(sizeof(struct LINK));
		(*g)->node = nodes+k;
		(*g)->length = len;
		(*g)->crash = cr;
		(*g)->more = NULL;
		
		for (g = &nodes[k].link; *g != NULL; g = &(*g)->more)
			;
		*g = malloc(sizeof(struct LINK));
		(*g)->node = nodes+j;
		(*g)->length = len;
		(*g)->crash = cr;
		(*g)->more = NULL;
	}
	
	nd_num = n; /* real number of nodes */

	return 0;

}

/* calculate distance matrix for crashes, using the linked list */

int dist_matrix()
{
	int i, j;

	struct LINK *s;
	struct CRASH *r;

	/* allocate memory for cr_dist */
	cr_dist = malloc(sizeof(double *) * cr_num);
	*cr_dist = malloc(sizeof(double) * cr_num * cr_num);

	/* allocate memory for cr_euc */
	cr_euc = malloc(sizeof(double *) * cr_num);
	*cr_euc = malloc(sizeof(double) * cr_num * cr_num);

	/* allocate memory for cr_det */
	cr_det = malloc(sizeof(double *) * cr_num);
	*cr_det = malloc(sizeof(double) * cr_num * cr_num);

	for (i=1; i < cr_num; i++)
	{	cr_dist[i] = cr_dist[i-1] + cr_num;
		cr_euc[i] = cr_euc[i-1] + cr_num;
		cr_det[i] = cr_det[i-1] + cr_num;
	}

//	printf ("-1 memory allocated..\n");

	/* fill cr_dist with very big values */
	for(i = 0; i < cr_num; i++)
		for(j = 0; j < cr_num; j++)
			cr_dist[i][j] = tot_len ;

//	printf ("-2 default values placed..\n");

	/* for every crash in the array of crashes */
	for (i = 0; i < cr_num; i++)
	{

//		printf ("\n- %d pt: ", i);

		/* calculate planar distances */
		for (j = 0; j < cr_num; j++)
			cr_euc[i][j] = dist(crashes[i].pt, crashes[j].pt);
		
//		printf ("pl, ");

		/* clean all the nodes.tot_dist */
		for(j = 0; j < nd_num; j++)
			nodes[j].tot_dist = 0.0;

//		printf ("cl, ");


		/* find the link where the crash is */
		
	

		if (crashes[i].adj_road == NULL)
			continue;
		else 
			for (s = crashes[i].pf->link; s->node != crashes[i].pl; s = s->more)
			;


		/* check the crashes on the same road */
		for (r = s->crash; r != NULL; r = r->more)
		{	cr_dist[i][r - crashes] = fabs(r->pf_dist - crashes[i].pf_dist);
			cr_dist[r - crashes][i] = fabs(r->pf_dist - crashes[i].pf_dist);
		}


		/* check the crashes on the other roads */
		crashes[i].pf->tot_dist = crashes[i].pf_dist;
		crashes[i].pl->tot_dist = crashes[i].pl_dist;
		go(crashes[i].pf, i);


		go(crashes[i].pl, i);


	}

	/* calculate detour index for each pair of points */
	for (i = 0; i < cr_num; i++)
		for (j = 0; j < cr_num; j++)
		{
			cr_det[i][j] = (double)(cr_euc[i][j]/cr_dist[i][j]);
			if (i == j) cr_det[i][j] = 1;
			if (cr_det[i][j] > 1) cr_det[i][j] = 1;
		}

	return 0;

}


/* calculate the total amounts (and indices) */

int result()
{
	int i, j, n, m, m1, dn;
	double dh, di=0, n1=0;

	/* define the sizes of the arrays */

	w = malloc(sizeof(int) * (int)((max_dist + d - st) / d + 1)); /* check, sometimes unexpected eresults */
	k = malloc(sizeof(int) * (int)((max_dist + d - st) / d + 1));
	b = malloc(sizeof(double) * (int)((max_dist + d - st) / d + 1));
	kb = malloc(sizeof(double) * (int)((max_dist + d - st) / d + 1));
	p = malloc(sizeof(int) * (int)((max_dist + d - st) / d + 1));
	p1 = malloc(sizeof(int) * (int)((max_dist + d - st) / d + 1));
	pk = malloc(sizeof(double) * (int)((max_dist + d - st) / d + 1));
	dt = malloc(sizeof(double) * (int)((max_dist + d - st) / d + 1));

	det = malloc(sizeof(double) * cr_num);

	
	/* print average detour for each crash */
	fp1 = fopen("detour.txt", "w");
	m = 0;
	for(i = 0; i < cr_num; i++)
	{	
		if (crashes[i].adj_road == NULL)
			continue;
		else 
		{
		det[i] = 0;
		n = 0;
		for(j = 0; j < cr_num; j++)
			if(cr_dist[i][j] != tot_len && cr_det[i][j] > 0) /* work on accuracy!! */
			{	
				det[i] += cr_det[i][j];
				n += 1;
			}

//		fprintf(fp1, "\n");	

		/* put in a file the values of average detour indices */
		fprintf(fp1, "%10.3f\n", (det[i] / n));

		di+=(det[i]/n);
		}
	}
	di = di / (cr_num - cr_null);
	fprintf(fp1, "Average detour index = %10.3f\n", di);


	/* create the output file. TODO: check if it is alrready created!! */
	
	fp = fopen("output.txt", "w");

	fprintf(fp, "Network threshold, m\tCount by network\tCount per bin by network\tNetwork K-function\tNetwork K-function per bin\tCount by plane(%5.3f of network threshold)\tCount by plane\tPlanar K-function\tDetour Index\n", di);

	for (h = st; h < max_dist + d; h+=d) /* for each threshold */
	{	a = (int)((h-st)/d);
		n = 0; /* counter */
		m = 0;
		m1 = 0;
		dh=0;
		dn=0;
		for(i = 0; i < cr_num; i++)
			for(j = 0; j < cr_num; j++)
			{	if (cr_dist[i][j] < h)
				{	n++;
					if(cr_det[i][j] > 0) 
					{	dh += cr_det[i][j];
						dn++;
					}
				}
				if (cr_euc[i][j] < (di * h))
					m++;
				if (cr_euc[i][j] < (h))
					m1++;

			}

		/* calculate values */
		w[a] = (int)(n / (cr_num - cr_null)); /* average amount of crashes within the network distance */
		k[a] = (int)((tot_len * w[a])/(cr_num - cr_null-1)); /* network K-function */
		b[a] = (double)((n - n1)/(cr_num - cr_null)); /* average amount of crashes per current bin */

		if (h==st)
			b[a]--; 

		kb[a] = (double)((tot_len*b[a])/(cr_num - cr_null -1)); /* network K-function per current bin */

		p[a] = (int)(m / cr_num); /* average amount of crashes within the corresponding Euclidean distance */
		p1[a] = (int)(m1 / cr_num); /* average amount of crashes within the same Euclidean distance */
		pk[a] = (double)((100000*p1[a])/cr_num); /* planar K-function */

		dt[a] = (double)(dh/dn); /* detour within this network distance */
		n1 = n;

		/* put in a file the value of threshold, tot. amount of crashes, crashes per bin and K-function */
		fprintf(fp, "%10.3f\t%10d\t%10.3f\t%10d\t%10.3f\t%10d\t%10d\t%10.3f\t%10.3f\n", h, w[a], b[a], k[a], kb[a], p[a], p1[a], pk[a], dt[a]);

	}

	fclose (fp);

	return 0;
}


/* read the point coordinates from the string of .txt file */

int get_coord(FILE *fp, double *x, double *y)
{
	int i=0;
	int k=1;
	double pow=0.0;
	
	*x = 0.0;
	*y = 0.0; 

	for (a = fgetc(fp); a != ')' ; a=fgetc(fp)) 
	{	/* read the coordinates */
		if isdigit(a) 
			if (k == 1)
			{	*x = 10.0 * *x + (a - '0');
				pow *= 10;
			}else if (k == 2)
				{*y = 10.0 * *y + (a - '0');
				pow *= 10;
				}
		if (a == '.')
			pow = 1.0;
		if (a == ',')
		{	if (k == 1)
			{	*x = *x / pow; 
				k++;
			} else
			if (k == 2)	
			{	*y = *y / pow;
				break; 
			}
		}
	}
	return 0;
}


/* create a new segment in the road from coordinates in .txt file */
/* if lines are broken - this part is not used */

struct ROAD *add_segment(FILE *fp, struct ROAD *segment) /* asks the pointer on the root segment */
{	
	/* allocate memory for new segment */
	segment->link = malloc(sizeof(struct ROAD));
	
	segment->link->pf = segment->pl;
	get_coord(fp, &segment->link->pl.x , &segment->link->pl.y );
	segment->link->link = NULL;
	segment->link->crashes = NULL;
	
	while ( (a=fgetc(fp)) != '\n' ) 
		/* skip all white and unnecessary symbols until the end of line */;

	for (a = fgetc(fp); a == ' ' || a == '	' ; a=fgetc(fp))
		;
	
	if (a == '(')
		segment->link->link = add_segment(fp, segment->link);
	
	return segment->link;
}



/* check if the crash is within any of the segments of the polyline */

int crash_segment(int i, int j, struct ROAD *segm)
{
	struct CRASH *q, *p;
	struct ROAD *z;

	/* check segment if the crash is within its extent+70m */
	if (crashes[j].pt.x <= max (segm->pf.x, segm->pl.x) + 70.0 && crashes[j].pt.x
	 >= min (segm->pf.x, segm->pl.x) - 70.0 && crashes[j].pt.y <= max (segm->pf.y,
	 segm->pl.y) + 70.0 && crashes[j].pt.y >= min (segm->pf.y, segm->pl.y) - 70.0)
	{
		/* if the distance beween the crash and the segment is less than saved for the crash */
		if (fabs(pt_dist(crashes[j].pt, *segm))< crashes[j].r_dist)
		{	
			/* 1) erase the record about the crash in the previous adjasent road */
			
			for (z = crashes[j].adj_road; z != NULL; z = z->link)
				if (z->crashes == crashes + j)
				{	p = z->crashes->more;
					z->crashes->more = NULL;
					z->crashes = p;
					break;
				}else if (z->crashes != NULL)
				{	for (q = z->crashes; q->more != crashes + j && q->more != NULL; q = q->more)
						;
					if (q->more == crashes + j)
					{	p = q->more->more;
						q->more->more = NULL;
						q->more = p;
						break;
					}
				}

			/* 2) set the crash to the road */
			if (segm->crashes == NULL)
			{
				segm->crashes = crashes + j;
			}
			else 
			{
				for (q = segm->crashes; q->more != NULL; q = q->more) 
				q->more = crashes + j;

			}

			/* 3) rewrite the distance and set the new road to the crash */
			crashes[j].r_dist = fabs(pt_dist(crashes[j].pt, *segm));
			crashes[j].adj_road = roads + i; 
		}
	}
	return 0;
}


/* distance between a point and a line */ 
double pt_dist(struct POINT obj, struct ROAD segm)
{
	double k = (segm.pf.y - segm.pl.y)/(segm.pf.x - segm.pl.x);
	double b = segm.pf.y - k*segm.pf.x;
	double K = -1/k;
	double B = obj.y - K*obj.x;
	double X = (B-b)/(k-K);
	double Y = k*X +b;


	/* if the point is within the segment - perpendicular */

	if (obj.x <= (max (segm.pf.x, segm.pl.x) + 5.0) && obj.x >= (min (segm.pf.x, segm.pl.x) - 5.0) && obj.y <= (max (segm.pf.y,
	 segm.pl.y) + 5.0) && obj.y >= (min (segm.pf.y, segm.pl.y) -5.0))
	{
		return fabs(((segm.pf.y-segm.pl.y)*obj.x + (segm.pl.x-segm.pf.x)*obj.y + (segm.pf.x*segm.pl.y - segm.pl.x*segm.pf.y)) / dist(segm.pf, segm.pl)); 
	}
	else /* else - distance to the closest end */
		return min (dist(obj, segm.pf), dist (obj, segm.pl));

}


/* distance between two points - Pifagor theorem */

double dist(struct POINT p1, struct POINT p2)
{
	return sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
}


/* go to the adjasent nodes, calculate the tot_dist, if the recorded is bigger */
/* also check for crashes */
/* algorythm now is closer to the depth-search, maybe it is needed to change it */

int go(struct NODE *n, int i)
{
	struct LINK *k;
	struct CRASH *c;

	/* go to all adjasent nodes (links) */
	for (k = n->link; k != NULL; k = k->more)
	{
		/* check, to go further or not */
		if (k->node->tot_dist > n->tot_dist + k->length || k->node->tot_dist == 0.0)
		{	k->node->tot_dist = n->tot_dist + k->length;
			if (k->node->tot_dist < max_dist)
				go(k->node, i);
		}

		/* look for a crash if the distance is less -> write */
		for (c = k->crash; c != NULL; c = c->more)
		{	if (c->pf == n)
				if (cr_dist[i][c - crashes] > n->tot_dist + c->pf_dist)
				{	cr_dist[i][c - crashes] = n->tot_dist + c->pf_dist;
					cr_dist[c - crashes][i] = n->tot_dist + c->pf_dist;
				}
			else 
				if (c->pl == n)
					if (cr_dist[i][c - crashes] > n->tot_dist + c->pl_dist)
					{	cr_dist[i][c - crashes] = n->tot_dist + c->pl_dist;
						cr_dist[c - crashes][i] = n->tot_dist + c->pl_dist;
					}
		}		
	}
	return 0;
}
