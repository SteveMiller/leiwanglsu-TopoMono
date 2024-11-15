//#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
//#include <windows.h>
//#include <tchar.h>
#include <stdlib.h>
//#include <bool.h>
#include "DrainageNetwork.h"

// I/O files defined
FILE *gfpSDP_Vtxs;					// The text file containing the ordered vertices of the Surface Drainage Path lines - gAllVtxs
FILE *gfpDEMAsPoints;				// The LiDAR DEM stored as xyz (*.tif file exported as x , y, z)
FILE *gfpSDP_VtxsLogFile;			// The log-file for processing results
FILE *gfpOrderedPlys;				// Final properly-ordered gAllPolylines	
FILE *gfpOrderedVtxs;				// Final properly-ordered gAllVtxs
FILE *gfpQuededFeatures;			// A list of erroneous polylines and nodes that require review and possible edit
FILE *gfpDrainFeaturesNetwork;		// The network table for polylines to nodes

// Note: gfpOrdered features plus gfpDrainFeaturesNetwork = Surface Drainage Network (SDN)

// GLOBALS - heavy reliance on global variables follows - dangerous but it clears up the function stacks

// A whole line of input from the tab-delimited text file; used by fgets()
char gLine[1024];

#define TRUE 1
#define FALSE 0

// Global DEBUG variables used to set and examine contents of multi-levelled pointers
struct VTX *gLookVTX;
struct POE *gLookPOE; 
struct PLY *gLookPLY;
struct NOD *gLookNOD;

// Limits on the DEM
double gDEM[1001][1001];
double gDEMxOrig = 308000.5;
double gDEMyOrig = 4408999.5;
double gDEMxySpace = 1.0; 

// Parameters for the PEM - redefined as a VTX_LST. Note the shift in origin.
double gPEMxOrig = 308000.0;
double gPEMyOrig = 4409000.0;
double gPEMxySpace = 1.0;
	
int gNumRows =1000, gNumCols =1000, gNumPosts;

struct SEN_VTX_LST *gPEM;

int gNumFields;

struct VTX *gAllVtxs;
struct VTX *gCurrVTX;
int gNumVTXs = 0;

struct PLY *gAllPolylines;
struct PLY *gCurrPLY;
int gNumPLYs = 0;

struct NOD *gNetworkTable;
struct NOD *gCurrNOD;
int gNumNODs = 0;

struct PLY *gListOfPLYErrors[20000];
int gNumPLYErrors = 0;

int gNumNetworkMismatches = 0;
int gNumReversed = 0;

double gThresholdDiff = 0.0;
struct PLY *gPLYsBelowThreshold[5000];
int gNumThreshPLYs = 0;

//******************************************************************************
// Functions defined below
//******************************************************************************
	
double BiLinearEstimateOfZ(double x, double y)
{
	int i, j;
	double Row, Col;
	double zNW, zNE, zSE, zSW;
	double A, B, C, D;
	double wtNW, wtNE, wtSE, wtSW;
	double z;

	// Find the floating point estimates of row, column at x, y
	Row = (gDEMyOrig - y) / gDEMxySpace;
	Col = (x - gDEMxOrig) / gDEMxySpace;
	
	// Convert floating point row, column to integer row, column (i,j)
	i = (int)Row;
	j = (int)Col;
	
	// Fetch the z-values at the bounding corners of the point's neighborhood
	zNW = gDEM[i][j];
	zNE = gDEM[i][j+1];
	zSE = gDEM[i+1][j+1];
	zSW = gDEM[i+1][j];
	
	// Simplify the calculations for weights
	A = 1.0 - (Row - (int)Row);
	B = 1.0 - (Col - (int)Col);
	C = Row - (int)Row;
	D = Col - (int)Col;

	// Calculate the weights for z at x, y
	wtNW = A * B;
	wtNE = A * D;
	wtSE = C * D;
	wtSW = C * B;

	// Apply the weights to derive a bilinear estimate of z at x, y
	z = wtNW * zNW + wtNE * zNE + wtSE * zSE + wtSW * zSW;
	
	return(z);
}	
	
int LoadDEMPosts()
{
	int i,j;
	int NumFields;
	double x, y;
	char xIn[25], yIn[25],zIn[25];
	
	// Skip past the header
	fgets(gLine, 1024, gfpDEMAsPoints);
	
	for(i =0;i <1000; i++)
	{
		gNumRows = gNumRows + 1;
		for(j=0; j <1000; j++)
		{
			gNumCols = gNumCols + 1;
			fgets(gLine, 1024, gfpDEMAsPoints);
			NumFields = sscanf(gLine, "%s\t%s\t%s", &xIn, &yIn, &zIn);
			// Ignore x and y for present, only record the z for the DEM
			gDEM[i][j] = atof(zIn);
		}
	}
	return (i*j);
}

int PopulateVertex(int *Network)
{
	// gfpSDP, gfpOWB, and gfpDLD file contents look like this:
	//  X			Y		Z	ID	Net	IDX		Dst		Ang
	//308996.5	4408188.5	0	3	1	0	0			135
	//308997.5	4408187.5	0	3	1	1	1.414213562	135
	//308998.5	4408186.5	0	3	1	2	2.828427125	135
	//308999.5	4408185.5	0	3	1	3	4.242640687	135
	//308997.5	4408132.5	0	4	2	0	0			45
	//308998.5	4408133.5	0	4	2	1	1.414213562	45
	//308999.5	4408134.5	0	4	2	2	2.828427125	45

	int Flag = 0;
	
	char xIn[15], \
		 yIn[15], \
		 zIn[15], \
		 idIn[15], \
		 networkIn[15],\
		 vertex_indexIn[15], \
		 distanceIn[15], \
		 angleIn[15];

	// While we still have a new record, load attributes of the VTX...
	// Read every time, but only store the Network data once on a PLY, not the VTXs
	// Ignore the indexIn - recompute vertex index values
	// Ignore the angle for now
	if(fgets(gLine, 1024, gfpSDP_Vtxs) != NULL)
	{
		gNumFields = sscanf(gLine, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
			&xIn, \
			&yIn, \
			&zIn, \
			&idIn, \
			&networkIn, \
			&vertex_indexIn, \
			&distanceIn, \
			&angleIn); 

		//Check to see if this was a valid sscanf
		if(gNumFields == 8)
		{	
			// Assign the vertex using gAllVtxs + gNumVTXs
			gCurrVTX = (struct VTX *)(gAllVtxs + gNumVTXs);

			// Populate the VTX data
			gCurrVTX->x = atof(xIn);
			gCurrVTX->y = atof(yIn);
		
			// In the export of the SDP polyline data from shapefiles, z-values were dropped.  Now we need to interpolate z from the DEM and x, y of the VTX.
			// Note that the bi-linear estimate used here will replicate z precisely where x and y are integrals of post spacing as well as handle
			// "floating point" x and y where data has been manually captured and vertices do not lie precisely on post x and y.
			gCurrVTX->z = BiLinearEstimateOfZ(gCurrVTX->x, gCurrVTX->y);

			// Looks like the Network is an identifier for a stream network that may prove useful.  Values range from 1 - 189 meaning that there are 189
			// unique networks.  Most are very small and only 20 or so are visible in color-coded displays at scale and using Network as a display variable.
			// Return this value so it can be associated with the CurrPLY and not the individual vertices.
			*Network = atoi(networkIn);

			// Store cumulative distance upstream in meters initially - convert to normalized values (0.0 to 1.0) for previous PLY back in LoadSDPVtxs()
			// when a new VtxDist = 0 is encountered
			gCurrVTX->NormDist = atof(distanceIn);
		
			
			if(gCurrVTX->NormDist != 0.0)
			{
				Flag = 0;	// Succesful read of next consecutive VTX
			}
			else
			{
				Flag = 1;	// Encountered the beginning of a new PLY
			}
		}
		// Failed to read 8 fields - set Flag to 2
		else			
		{
			printf(" Failed to successfully sscanf 8 fields for Vertex %d", gNumVTXs);
			gCurrVTX = NULL;
			Flag = 2;	// Failed to scan 8 fields
		}
		
	}
	// If here, we have encountered an EOF on fgets, set Flag to 3
	else
	{
		printf(" Reached EOF on attempted read for Vertex %d", gNumVTXs);
		gCurrVTX = NULL;
		Flag = 3;
	}
	return(Flag);
}	
int LoadSDPVtxs()
{
	// gfpSDP, gfpOWB, and gfpDLD file contents look like this:
	
	//  X			Y		Z	ID	Net	IDX		Dst		Ang
	//308996.5	4408188.5	0	3	1	0	0			135
	//308997.5	4408187.5	0	3	1	1	1.414213562	135
	//308998.5	4408186.5	0	3	1	2	2.828427125	135
	//308999.5	4408185.5	0	3	1	3	4.242640687	135
	//308997.5	4408132.5	0	4	2	0	0			45
	//308998.5	4408133.5	0	4	2	1	1.414213562	45
	//308999.5	4408134.5	0	4	2	2	2.828427125	45
	
	int Flag;
	int Network, VtxOwnerID, VtxIndex;
	
	struct VTX *PrevPLYVTX;
	struct VTX *PrevVtx;
	struct PLY *PrevPLY;
	 

	// Skip past the header line for the vertex data
	fgets(gLine, 1024, gfpSDP_Vtxs);
	
	//*******************************************************************************************************************
	// In a single pass, read and sscanf for all input vertices creating a doubly-linked list of vertices which is then 
	// divided into polylines by severing next and previous links between the vertices as we go.  Polyline count and content
	// are updated when the current vertex has a NormDist == 0.0 indicating the start of a new PLY.
	//*******************************************************************************************************************/
	
	// Loop and process until we reach EOF
	while((Flag = PopulateVertex(&Network)) < 2)	// Either another consecutive VTX in the current PLY or a new PLY altogether
	{
		// Assign the PLY Data if NormDist was 0.0 - first VTX in a new PLY
		if(Flag == 1)
		{
		
			// Assign a current PLY
			gCurrPLY = (struct PLY *)(gAllPolylines + gNumPLYs);
			gCurrPLY->ID = gNumPLYs;	
			
			// Safe to assign a root vertex and a next VTX pointer
			gCurrPLY->RootVtx = gCurrVTX;
			
			// And safe to set pointer from this VTX to the next - as long as we stay in allocated memory bounds
			if(gNumVTXs < 249998)
			{
				gCurrVTX->NextVtx = (struct VTX *)(gAllVtxs + gNumVTXs + 1);
			}
			else
			{
				printf("\n Out of allocated memory limits for vertices. \n");
				Flag = 4;
				return(Flag);
			}
			
			// Safe to set previous VTX pointer
			gCurrVTX->PrevVtx = NULL;
			
			// For all but the first PLY...
			if(gNumPLYs > 0)
			{
				// Break the link from this, the first VTX in the new PLY, to the previous VTX next VTX
				PrevVtx = gCurrVTX - 1;
				PrevVtx->NextVtx = NULL;
			}
			
			// Update OwnerLine info for the VTX
			gCurrVTX->OwnerPLY = gCurrPLY;
	
			// Set the new PLY RootVtx
			gCurrPLY->RootVtx = gCurrVTX;
			
			// Look ahead a bit in processing and store the Network Index off the current Polyline structure. This saves space since
			// we do not have to record the redundant network value on each Vertex - now appropriately just recorded on the Polyline.
			gCurrPLY->NetworkIndex = Network;
			
			// Finish processing for the current PLY when we encounter a second new PLY (gNumPLYs > 0) and for all subsequent PLYs
			if(gNumPLYs > 0)
			{	
				// Set Term VTX for the previous PLY to previous VTX
				PrevPLY = (struct PLY *)(gAllPolylines + gNumPLYs - 1);
				PrevPLY->TermVtx = (struct VTX *)(gCurrVTX - 1);
				
				// Assign Previous PLY length from gCurrTX distance which is currently in meters
				PrevPLY->Length = PrevPLY->TermVtx->NormDist;
				

				gCurrPLY->RootVtx = (struct VTX *)gCurrVTX;								

				// Sever the link from last VTX in previous PLY to CurrPLY->RootVTX
				PrevPLY->TermVtx->NextVtx = NULL;
				
				// Now loop through all VTXs in the previous PLY computing NormDist
				PrevPLYVTX = (struct VTX *)(PrevPLY)->RootVtx;
				
				while (PrevPLYVTX != NULL)
				{
					// Use the actual distance stored in NormDist currently, divided by the LastVtx->NormDist, to render a NormDist in the range of 0.0 to 1.0
					PrevPLYVTX->NormDist = PrevPLYVTX->NormDist / (PrevPLY->Length);
					PrevPLYVTX = PrevPLYVTX->NextVtx;
				}					
			}	// End processing VTXs in the current PLY
			
			// Bump the offset index for a new PLY 
			gNumPLYs = gNumPLYs + 1;
			
		}	// End start of next Polyline
		// Flag = 0. This is an intermediate VTX.
		else
		{
			if(gNumVTXs < 249998)
			{
				if (gNumVTXs > 0)
				{
					// Link backwards and forward for all intermediate VTXs
					gCurrVTX->PrevVtx = (gAllVtxs + gNumVTXs - 1);
					gCurrVTX->NextVtx = (gAllVtxs + gNumVTXs + 1);
				}
				else
				{
					gCurrVTX->PrevVtx = NULL;
					gCurrVTX->NextVtx = (gAllVtxs + gNumVTXs + 1);
				}

				// Update OwnerLine info for the VTX
				gCurrVTX->OwnerPLY = gCurrPLY;	
			}
			else
			{
				printf("\n Out of allocated memory limits for vertices. \n");
				Flag = 4;
				return(Flag);	
			}
		}
		
		// Bump the VTX index
		gNumVTXs = gNumVTXs + 1;
		
	}	// EOF on data

	// We fell out of the read data loop with EOF true and now must process the current (last) PLY	

	// Reset current VTX to the last valid VTX in the PLY
	gCurrVTX = (struct VTX *)(gAllVtxs + gNumVTXs - 1);
	
	// Sever the last PLY's gCurrVTX Next
	gCurrVTX->NextVtx = NULL;
	
	// Set length of PLY to last cumulative distance
	gCurrPLY->Length = gCurrVTX->NormDist;
	
	// Set gCurrPLY->TermVtx
	gCurrPLY->TermVtx = gCurrVTX;
	
	// Loop through all vertices in the last polyline and re-compute NormDist
	PrevPLYVTX = (struct VTX *)(gCurrPLY)->RootVtx;
	
	while (PrevPLYVTX != NULL)
	{
		// Use the actual distance stored in NormDist currently, divided by the TermVtx->NormDist, to render a NormDist in the range of 0.0 to 1.0
		PrevPLYVTX->NormDist = PrevPLYVTX->NormDist / (gCurrPLY->Length);
		PrevPLYVTX = PrevPLYVTX->NextVtx;
	}
	
	printf("\n Loaded %d SDP vertices and %d Polylines \n", gNumVTXs, gNumPLYs);
}

void LoadOWBVtxs()
{
	// 100 LOC
}

void LoadDLDVtxs()
{
	// 200 LOC
}

void AddLinkToPLY_LST()
{
	
}

void CreateNODTable()
{
	int idx;
	
	gNetworkTable = (struct NOD *)calloc(100000, sizeof(struct NOD));
	
	for (idx = 0; idx < 100000; idx++)
	{
		(gNetworkTable + idx)->ID = idx;
		(gNetworkTable + idx)->DownstreamPolylines = NULL;
		(gNetworkTable + idx)->UpstreamPolylines = NULL;
	}
}

int CreatePEM()
{
	int i, j, index;
	
	gNumRows = 1000;
	gNumCols = 1000;
	
	gPEM = (struct SEN_VTX_LST *)calloc(gNumRows * gNumCols, sizeof(struct SEN_VTX_LST));
			
	for(i =0; i < gNumRows; i++)
	{	
		for(j=0; j < gNumCols; j++)
		{
			index = i * gNumCols + j;
			(gPEM + index)->HeadList = NULL;
			(gPEM + index)->TailList = NULL;
		}
	}
	return (i*j);
}

void EvaluateThresholdElevDiff()
{
	int i;
	int i_root, i_term, j_root, j_term;
	double Mean, Min, Max, Sum, Total, xDif, yDif, xyDif;
	double *Diff, *Slope, *SumSqd, *ZScore;
	struct PLY *Polyline;
	
	Max = -99999.0;
	Min = 99999.0;
	Diff = (double *) calloc(gNumPLYs, sizeof(double));
	SumSqd = (double *) calloc(gNumPLYs, sizeof(double));
	ZScore = (double *) calloc(gNumPLYs, sizeof(double));
	Slope = (double *) calloc(gNumPLYs, sizeof(double));
	
	i = 0;
	Sum = 0.0;
	
	while (i < gNumPLYs)
	{
		Polyline = (struct PLY *)(gAllPolylines + i);
	
		i_root = (gPEMyOrig - (Polyline->RootVtx)->y) / gPEMxySpace;
		j_root = (Polyline->RootVtx)->x - gPEMxOrig / gDEMxySpace;

		i_term = (gPEMyOrig - (Polyline->TermVtx)->y) / gPEMxySpace;
		j_term = (Polyline->TermVtx)->x - gPEMxOrig / gPEMxySpace;
	
		Diff[i] = fabs(gDEM[i_term][j_term] - gDEM[i_root][j_root]);
		
		if(Diff[i] < Min)
		{
			Min = Diff[i];
		}	
		
		if(Diff[i] > Max)
		{
			Max = Diff[i];
		}
		Sum = Sum + Diff[i];	
		i = i + 1;	
	}
	
	Mean = Sum / gNumPLYs;
	
	i = 0;
	Total = 0.0;
	
	while (i < gNumPLYs)
	{
		SumSqd[i] = pow(Diff[i] - Mean, 2);
		Total = Total + SumSqd[i];
		i = i + 1;
	}
	
	gThresholdDiff = sqrt(Total / (gNumPLYs - 1));
	
	i = 0;
	
	while (i < gNumPLYs)
	{
		ZScore[i] = (Diff[i] - Mean) / gThresholdDiff;
		xDif = (gAllPolylines + i)->TermVtx->x - (gAllPolylines + i)->RootVtx->x;
		yDif = (gAllPolylines + i)->TermVtx->y - (gAllPolylines + i)->RootVtx->y;
		xyDif = sqrt(xDif * xDif + yDif * yDif);
		Slope[i] = Diff[i] / xyDif;
		fprintf(gfpQuededFeatures,"%lf\t %lf \t %lf \n", Diff[i], xyDif,  Slope[i]);
		i = i + 1;
	}
	
	printf(" For %d gAllPolylines Mean = %lf Min = %lf Max = %lf and gThresholdDiff = %lf \n", gNumPLYs, Mean, Min, Max, gThresholdDiff); 
}

_Bool CompareZsOfAdjoiningPolylines(struct PLY *NextDownPLY, struct PLY *NextUpPLY)
{
	_Bool Match;
	
	struct VTX *DownVTXPtr;
	struct VTX *UpVTXPtr;
	
	DownVTXPtr = (struct VTX *)NextDownPLY->RootVtx;
	UpVTXPtr = (struct VTX *)NextUpPLY->TermVtx;
	
	// Is the downstream polyline's root vertex z less than or equal to the upstream polyline's term vertex z? If so, OK and Match is TRUE.
	if(DownVTXPtr->z <= UpVTXPtr->z)
	{
		Match = TRUE;
	}
	else
	{
		Match = FALSE;
	}
	return(Match);
}

int PiecewiseValidationOfGlobalNetworkConnectivity()
{
	_Bool Match;
	int idx, NumberMismatches;	
	struct NOD *CurrNODPtr;
	struct PLY_LST *DownLineNbrs;
	struct PLY_LST *UpLineNbrs;
	struct PLY *NextDownPLY;
	struct PLY *NextUpPLY;
	
	idx = 0;
	
	CurrNODPtr = (struct NOD *)(gNetworkTable + idx);
	
	// While we still have an unexamined node...
	while(idx < gNumNODs)
	{
		Match = TRUE;
		
		// A list of all polyline(s)downstream from the current node
		DownLineNbrs = (struct PLY_LST *)(CurrNODPtr->DownstreamPolylines);
		
		// The first polyline in a possible seqence of polylines downsream from the current node
		NextDownPLY = (struct PLY *)(DownLineNbrs->ThisPLY);
		
		// Do we have a downsteam polyline?
		while(NextDownPLY != NULL)
		{
			// A list of all polylines upstream from the current node
			UpLineNbrs = (struct PLY_LST *)(CurrNODPtr->UpstreamPolylines);			
			
			// The first polyline in a possible seqence of polylines upsream from the current node
			NextUpPLY = (struct PLY *)(UpLineNbrs->ThisPLY);

			// Do we have an upstream polyline?
			while(NextUpPLY != NULL)
			{
				// Compare Root VTX z of NextDownPLY with Term VTX z of NextUpPLY - are they monotonic?
				Match = CompareZsOfAdjoiningPolylines(NextDownPLY, NextUpPLY);
				
				// If no Match, then the lines are not montonic - queue the polylines for edit
				if(Match == FALSE)
				{
					NumberMismatches = NumberMismatches + 1;
					
					fprintf(gfpQuededFeatures, "Mismatch %d for downstream ID = %d at x = %lf y = %lf and z = %lf and upstream ID = %d at x = %lf y = %lf and z = %lf \n",
					NumberMismatches,
					NextDownPLY->ID, NextDownPLY->RootVtx->x, NextDownPLY->RootVtx->y, NextDownPLY->RootVtx->z,
					NextUpPLY->ID, NextUpPLY->RootVtx->x, NextUpPLY->RootVtx->y, NextUpPLY->RootVtx->z,
					gNumNetworkMismatches = gNumNetworkMismatches + 1);
				}
				// Next upstream PLY if there is one
				UpLineNbrs = UpLineNbrs->Next;
				
				if(UpLineNbrs != NULL)
				{
					NextUpPLY = (struct PLY *)(UpLineNbrs->ThisPLY);
				}
				else
				{
					NextUpPLY = NULL;
				}
			}
			
			// Next downstream PLY if there is one
			DownLineNbrs = DownLineNbrs->Next;
			
			if(DownLineNbrs != NULL)
			{
				NextDownPLY = (struct PLY *)(DownLineNbrs->ThisPLY);
			}
			else
			{
				NextDownPLY = NULL;
			}
		}
		
		// Move to the next NOD
		idx = idx + 1;
		
		CurrNODPtr = (struct NOD *)(gNetworkTable + idx);
	}
}

int UpdateEditQueueForNetwork()
{
	// 60 LOC
}

_Bool OutputAll()
{
	_Bool Status;
	
	Status = TRUE;
	
	// 80 LOC
	return(Status);
}

void OutputOrderedFeature(int idx)
{
	int j;
	
	gCurrPLY = (struct PLY *)(gAllPolylines + idx);
		
	printf("\n In OutputOrderedFeature for Polyline %d \n", idx);
	
	fprintf(gfpOrderedPlys, "%d \t %d \t %lf \t %d \t %d \t %d \t %d \t %d \t %d \t %d \n", \
		gCurrPLY->ID, gCurrPLY->NetworkIndex, gCurrPLY->Length, (gCurrPLY->RootVtx-gAllVtxs), \
		(gCurrPLY->TermVtx-gAllVtxs), (gCurrPLY->DownstreamNode - gNetworkTable), \
		(gCurrPLY->UpstreamNode - gNetworkTable), gCurrPLY->OwnerSHL, \
		gCurrPLY->EntityLeft, gCurrPLY->EntityRight);
	
	gCurrVTX = gCurrPLY->RootVtx;
	
	j = -1;
	
	while(gCurrVTX != NULL)
	{
		j = j + 1;
		
		fprintf(gfpOrderedVtxs, "%d \t %d \t %lf \t %lf \t %lf \t %lf \n", \
			gCurrPLY->ID, j, gCurrVTX->x,gCurrVTX->y, gCurrVTX->z, gCurrVTX->NormDist);
		
		gCurrVTX = gCurrVTX->NextVtx;
	}
}

// Display the NOD to PLYs relations as they are established. Shows work in progress.
void OutputWorkingtNetworkTable(int idx)
{	
	struct PLY *DownPLY;
	struct PLY *UpPLY;
	struct PLY_LST *DownPLY_LST;
	struct PLY_LST *UpPLY_LST;
	struct PLY *ThisPLY;
	struct NOD *DownNOD;
	struct NOD *UpNOD;
	
	// Assign a local version of the Polyline
	ThisPLY = (struct PLY *)(gAllPolylines + idx);
	
	// 1. Check DOWN_NOD_TO_DOWN_PLYS. 
	// Checks for the Downstream Node for Downstream neighbor polylines
	// By now, some combination of Upstream and Downstream PLYs should have been established
	DownNOD = (struct NOD *)(ThisPLY->DownstreamNode);

	// No Downstream neighbors for Downstream NOD yet
	if(DownNOD->DownstreamPolylines != NULL)
	{
		DownPLY_LST = (DownNOD->DownstreamPolylines->HeadList);

	
		fprintf(gfpDrainFeaturesNetwork," Downstream Node ID %d \t to Down PLYs:", DownNOD->ID);
		
		while (DownPLY_LST != NULL)
		{
			DownPLY = DownPLY_LST->ThisPLY;
			fprintf(gfpDrainFeaturesNetwork,"\t %d", DownPLY->ID);
			DownPLY_LST = DownPLY_LST->Next;
		}
	}
	
	// No Downstream Polylines for the Downstream Node
	else
	{
		fprintf(gfpDrainFeaturesNetwork,"Downstream Node ID %d to Down PLYs \t NONE", DownNOD->ID);		
	}

	// 2. Check DOWN_NOD_TO_UP_PLYS. Checks for the Downstream Node and adjoining Upstream polylines
	if(DownNOD->UpstreamPolylines != NULL)
	{
		UpPLY_LST = (DownNOD->UpstreamPolylines->HeadList);

		fprintf(gfpDrainFeaturesNetwork," Downstream Node ID %d \t to Up LYs:", DownNOD->ID);
		
		while (UpPLY_LST != NULL)
		{
			UpPLY = UpPLY_LST->ThisPLY;
			fprintf(gfpDrainFeaturesNetwork,"\t %d", UpPLY->ID);
			UpPLY_LST = UpPLY_LST->Next;
		}
	}
	
	// No Upstream Polylines for the Upstream Node
	else
	{
		fprintf(gfpDrainFeaturesNetwork,"Upstream Node ID %d Up PLYs \t NONE", DownNOD->ID);		
	}
	
	// 3. Check UP_NODE_TO_UP_PLYS. Checks for the Upstream Node and adjoining Upstream polylines
	// Parallel check as for 1 above
	UpNOD = (struct NOD *)(ThisPLY->UpstreamNode);

	if(UpNOD->UpstreamPolylines != NULL)
	{
		UpPLY_LST = (UpNOD->UpstreamPolylines->HeadList);
	
		fprintf(gfpDrainFeaturesNetwork," Upstream Node ID %d \t Up PLYs:", UpNOD->ID);
		
		while (UpPLY_LST != NULL)
		{
			UpPLY = UpPLY_LST->ThisPLY;
			fprintf(gfpDrainFeaturesNetwork,"\t %d", UpPLY->ID);
			UpPLY_LST = UpPLY_LST->Next;
		}
	}
	
	// No Upstream Polylines for the Downstream Node
	else
	{
		fprintf(gfpDrainFeaturesNetwork,"Downstream Node ID %d Up PLYs \t NONE", DownPLY->ID);
	}
	
	// 4. Check UP_NODE_TO_DOWN_PLYS
	// Parallel check as for 2 above
	UpNOD = (struct NOD *)(ThisPLY->UpstreamNode);

	if(UpNOD->DownstreamPolylines != NULL)
	{
		DownPLY_LST = (UpNOD->DownstreamPolylines->HeadList);
	
		fprintf(gfpDrainFeaturesNetwork," Upstream Node ID %d \t Up PLYs:", UpNOD->ID);
		
		while (DownPLY_LST != NULL)
		{
			DownPLY = DownPLY_LST->ThisPLY;
			fprintf(gfpDrainFeaturesNetwork,"\t %d", DownPLY->ID);
			DownPLY_LST = DownPLY_LST->Next;
		}
	}
	
	// No Upstream Polylines for the Downstream Node
	else
	{
		fprintf(gfpDrainFeaturesNetwork,"Downstream Node ID %d Up PLYs \t NONE", DownPLY->ID);
	}
		
	// Print a carriage return to end the output line for all PLYs for both Up and Down NODs
	fprintf(gfpDrainFeaturesNetwork, " \n");
}

// Cycle through the vertices, swapping next and prev pointers for each vertex
void ReversePolyline()
{
	struct VTX *FirstVtx; 
	struct VTX *LastVtx;
	struct VTX *NextVtx;
	struct VTX *CurrVtx;
	struct VTX *SaveVtx;


	// Save positions of the endpoints to reverse Root and Term VTX pointers for the Polyline
	FirstVtx = (struct VTX *)gCurrPLY->TermVtx;
	LastVtx = (struct VTX *)gCurrPLY->RootVtx;
	
	// Set New Root VTX using current relations
	CurrVtx = FirstVtx;
	NextVtx = CurrVtx->PrevVtx;
	while (CurrVtx != NULL)
	{
		// Flip the relations for Curr VTX
		SaveVtx = CurrVtx->NextVtx;
		CurrVtx->NextVtx = CurrVtx->PrevVtx;
		CurrVtx->PrevVtx = SaveVtx;

		// Invert the normalized distance along the PLY for this VTX
		CurrVtx->NormDist = 1.0 - CurrVtx->NormDist;
		
		// Use the newly assigned Next to go "backwards" through the Polyline to establish the new direction for the PLY
		CurrVtx = CurrVtx->NextVtx;
	}

	gCurrPLY->RootVtx = FirstVtx;
	gCurrPLY->RootVtx->PrevVtx = NULL;
	gCurrPLY->TermVtx = LastVtx;
	gCurrPLY->TermVtx->NextVtx = NULL;
}

// MapNODs performs the following:
// 1. Computes i_root, j_root and i_term and j_term cell row and column in the PEM from Polyline endpoints x and y
// 2. Performs a quick assessment of Polyline orientation.  By convention, PLY VTXs should be ordered low to high in z and z(i+1) - z(i)  > 0
// 3. If z(i+1) - z(i) < 0.0, then call ReversePolyline()and reverse (re-compute) i_root, j_root and i_term and j_term
// 4. If endpoints ABS(z(i+1) - z(i)) < gThresholdDiff, the flow is considered indeterminant or ambiguous.  Put *PLY on an edit queue.
// 5. Assign PEM VTX_LSTs to Head and Tail LSTs
// 6. Assign PLYs to Upstream and Downstream LSTs for NODs
  
int MapNODs(int idx)
{
	int i_root, i_term, j_root, j_term;
	double zDif;
	struct PLY *ComparePLY;
	
	gCurrPLY = (struct PLY *)(gAllPolylines + idx);
	
	// Initial i_root, j_root, i_term and j_term computed using the original order of the PLY VTXs
	i_root = (int)((gPEMyOrig - gCurrPLY->RootVtx->y) / gPEMxySpace);
	j_root = (int)((gCurrPLY->RootVtx->x - gPEMxOrig) / gPEMxySpace);

	i_term = (int)((gPEMyOrig - gCurrPLY->TermVtx->y) / gPEMxySpace);
	j_term = (int)((gCurrPLY->TermVtx->x - gPEMxOrig) / gPEMxySpace);
	
	// A quick check for monotonicity and drain ambiguity by comparing zDif to difference statistics from endpoints of all Polylines
	// * NOTE - we may prefer to calculate zDif as Polyline->TermVtx->z - Polyline->RootVtx->z?!  Especially with Manual Data Capture option
	zDif = gDEM[i_term][j_term] - gDEM[i_root][j_root];
	
	// Threshold defaults to 0.05 * gThresholdDiff computed in EvaluateThresholdElevDiff()
	if ((zDif < 0.0) && (fabs(zDif) > (0.05 * gThresholdDiff)))
	{		
		// Anticipated case: well-defined drainage, but inverted high to low endpoints - reverse the PLY
		ReversePolyline();
		gNumReversed = gNumReversed + 1;
		
		// After the swap of RootVtx and TermVtx, re-compute the i_root, j_root, i_term, j_term values  
		i_root = (int)((gPEMyOrig - gCurrPLY->RootVtx->y) / gPEMxySpace);
		j_root = (int)((gCurrPLY->RootVtx->x - gPEMxOrig) / gPEMxySpace);

		i_term = (int)((gPEMyOrig - gCurrPLY->TermVtx->y) / gPEMxySpace);
		j_term = (int)((gCurrPLY->TermVtx->x - gDEMxOrig) / gPEMxySpace);
	}
	// Trap the ambiguous flow case using 0.05* gThresholdDiff - a small percentage of the standard deviation of all Polyline endpoint differences
	else if(fabs(zDif) <= (0.05 * gThresholdDiff))
	{
		// Add the PLY pointer to the array of ambiguous SDP polylines
		if(gNumThreshPLYs < 5000)
		{
			gPLYsBelowThreshold[gNumThreshPLYs] = gCurrPLY;
			gNumThreshPLYs = gNumThreshPLYs + 1;			
		}
		else
		{
			printf("Limit of 5000 Polylines below Threshold reached. \n STOP \n");
			return(1);
		}
	}
	else
	{
		// Do nothing - the zDif > (0.05 * gThresholdDiff) and this is a naturally well-ordered drain
	}

	// Debug writes to check output data.
	OutputOrderedFeature(idx);
	
	//*******************************************************************************************************************
	// Now, link the PEM, PLY and NODs for both root and term NODs of this Polyline!  Good luck, you're going to need it.
	//*******************************************************************************************************************	
	
	// First, the Root (Downstream) NOD "with respect to"(wrt)the current PLY...

	// There is a variable-length list of PLYs (a series of linked PLY_LST elements) associated with a NOD
	// HeadList and TailList that respectively point to the begining and end of the list
	// If the HeadList associated with gPEM at cell i_root, j_root is NULL, there are no currently assigned PEM entities.  
	// We need to 1) create initial gPEM entities for VTXs, 2)create a new upstream NOD, and 3) link the current PLY to the new NOD. 
	if((gPEM + i_root * gNumCols + j_root)->HeadList == NULL)
	{
		// Remember: each gPEM cell was allocated as a SEN_VTX_LST with Head and Tail pointers to the VTXs associated 
		// with the gPEM Cell.  Allocate a new HeadList and TaiList for the gPEM - each a pointer to the start and end VTX_LST 
		// associated with this gPEM element
		(gPEM + i_root * gNumCols + j_root)->HeadList = (struct VTX_LST *) calloc(1, sizeof(struct VTX_LST));
		(gPEM + i_root * gNumCols + j_root)->TailList = (struct VTX_LST *) calloc(1, sizeof(struct VTX_LST));	
		
		// Assign the Current PLY's RootVtx to both the HeadList and TailList's VTX
		(gPEM + i_root * gNumCols + j_root)->HeadList->ThisVTX = gCurrPLY->RootVtx;
		
		// Both Head and Tail point to the same VTX intially
		(gPEM + i_root * gNumCols + j_root)->TailList->ThisVTX = gCurrPLY->RootVtx;
		
		// For any first gPEM[i][j] HeadList entry, the previous VTX_LST pointer should be NULL
		// Again, setting the Tail of the List to the same value as the Head initially
		(gPEM + i_root * gNumCols + j_root)->HeadList->Prev = NULL;
		
		(gPEM + i_root * gNumCols + j_root)->TailList->Prev = NULL;
		
		// Set the Next PEM VTX LST to NULL as well - we will expand the list by expanding the Tail 
		// of the List as needed when this new NOD is later matched as seen in the "else" section below
		(gPEM + i_root * gNumCols + j_root)->HeadList->Next = NULL;

		(gPEM + i_root * gNumCols + j_root)->TailList->Next = NULL;
		
		// Need to assign a new NOD from the allocated space for NODs
		gCurrNOD = (struct NOD *)(gNetworkTable + gNumNODs);
		
		// Assign a VTX to the NOD from current Polyline
		gCurrNOD->ThisVTX= gCurrPLY->RootVtx;
		
		// Bump the offset index to gNetworkTable for next use
		gNumNODs = gNumNODs + 1;
			
		// Wrt to the current PLY, the current NOD is the downstream NOD
		// Set Polyline Downstream Node using the next available NOD structure from the Network Table
		gCurrPLY->DownstreamNode = gCurrNOD;

		// Wrt the current NOD, the current PLY is UPSTREAM of the NOD!
		// There may be several PLYs upstream (and/or downstream) of a NOD, therefore we have another 
		// variable-length list, this time for PLYs
		// Unlike the gPEM array of cells which were allocated upfront, here we create the Upstream SEN_PLY_LST for this NOD
		gCurrNOD->UpstreamPolylines = (struct SEN_PLY_LST *)calloc(1, sizeof(struct SEN_PLY_LST));
		
		// Allocate space for the Head and Tail PLY LSTs
		gCurrNOD->UpstreamPolylines->HeadList = (struct PLY_LST *)calloc(1, sizeof(struct PLY_LST));
		gCurrNOD->UpstreamPolylines->TailList = (struct PLY_LST *)calloc(1, sizeof(struct PLY_LST));	
		
		// Set Prev	pointers for the Head and Tail LSTs to NULL
		// Initially Head and Tail point to the single entry in the LST
		gCurrNOD->UpstreamPolylines->HeadList->Prev = NULL;
		gCurrNOD->UpstreamPolylines->TailList->Prev = NULL;
		
		// Set Next pointers for the Head and Tail LSTs to NULL
		// Initially Head and Tail point to the single entry in the LST
		gCurrNOD->UpstreamPolylines->HeadList->Next = NULL;
		gCurrNOD->UpstreamPolylines->TailList->Next = NULL;
		
		// Add the current Polyline to the Head and Tail PLY_LSTs for UpstreamPolylines 
		gCurrNOD->UpstreamPolylines->HeadList->ThisPLY = (struct PLY *)gCurrPLY;
		gCurrNOD->UpstreamPolylines->TailList->ThisPLY = (struct PLY *)gCurrPLY;
	}
	// Else, the HeadList is not NULL and a NOD already exists for this PEM[i][j]. The question becomes, 
	// how can we efficiently find this NOD?
	else
	{
		// A NOD already exists. Find the NOD by using a previously assigned entry and it's associated 
		// This VTX and Owner PLY. Use the TailList pointer to quickly find Compare PLY and append it to a 
		// new entry for Tail List.
		ComparePLY = (gPEM + i_root * gNumCols + j_root)->TailList->ThisVTX->OwnerPLY;

		// Confirm the NOD by comparing x and y of the current Polyline with those of the newly found ComparePLY
		// An exact check in xy works because coordinates in the SDP data are aligned with the grid post xy's precisely.
		if(gCurrPLY->RootVtx->x == ComparePLY->RootVtx->x && gCurrPLY->RootVtx->y == ComparePLY->RootVtx->y)
		{
			// Polylines match at their Downstream node; for each PLY - 1) assign PLY Downstream NOD and
			// 2)add current PLY to the list of UPstream PLYs for the NOD
			gCurrNOD = (struct NOD *)ComparePLY->DownstreamNode;
			
			// Wrt to the curent PLY, this is the Downstream NOD
			gCurrPLY->DownstreamNode = gCurrNOD;
			
			// Add the Curr PLY to the NOD's UPSTREAM(!) PLY_LST
			gCurrNOD->ThisVTX = gCurrPLY->RootVtx;
			
			// If we are just adding the second VTX to the PEM cell, special handling is needed for prev and next for the VTX_LSTs
			if((gPEM + i_root * gNumCols + j_root)->TailList->ThisVTX == (gPEM + i_root * gNumCols + j_root)->HeadList->ThisVTX)
			{
				// Initial Tail List update - re-use the existing Tail List
				(gPEM + i_root * gNumCols + j_root)->HeadList->Next = (gPEM + i_root * gNumCols + j_root)->TailList;
				(gPEM + i_root * gNumCols + j_root)->TailList->Prev = (gPEM + i_root * gNumCols + j_root)->HeadList;

				// Assign the current PLY Root VTX to the Tail List This VTX
				(gPEM + i_root * gNumCols + j_root)->TailList->ThisVTX = gCurrPLY->RootVtx;
			}
			// Else we are adding a third, fourth, fifth VTX...always extend from the Tail of the List
			else
			{
				// Update gPEM root TailList (a pointer to a VTX_LST structure)
				(gPEM + i_root * gNumCols + j_root)->TailList->Next = (struct VTX_LST *)calloc(1, sizeof(struct VTX_LST));	
				(gPEM + i_root * gNumCols + j_root)->TailList->Next->Prev = (gPEM + i_root * gNumCols + j_root)->TailList;
				
				// Update to a new List Tail pointer
				(gPEM + i_root * gNumCols + j_root)->TailList = (gPEM + i_root * gNumCols + j_root)->TailList->Next;

				// calloc() NULLs everything implicitly, but be deliberately explicit and set Next to NULL 
				(gPEM + i_root * gNumCols + j_root)->TailList->Next = NULL;				

				// Assign the current PLY Root VTX to the Tail List This VTX
				(gPEM + i_root * gNumCols + j_root)->TailList->ThisVTX = gCurrPLY->RootVtx;
			}			
		}
		// Does Current PLY's Root VTX match at the terminal node of Compare PLY? 
		// If so, we have matched one upstream and one downstream PLY at the NOD 
		else if(gCurrPLY->RootVtx->x == ComparePLY->TermVtx->x && gCurrPLY->RootVtx->y == ComparePLY->TermVtx->y)
		{
			// Matches at an upstream node for Compare PLY
			gCurrNOD = (struct NOD *)ComparePLY->UpstreamNode;

			// Note - while the Curr PLY Root VTX matched the Compare PLY Term VTX and it's Upstream NOD,
			// we need to assign the PLY to the NOD's Downstream LST. It may be necessary to draw a picture.
			// Extend the list of polylines and link the LST elements
			gCurrNOD->DownstreamPolylines->TailList->Next = (struct PLY_LST *)calloc(1, sizeof(struct PLY_LST));		
			gCurrNOD->DownstreamPolylines->TailList->Next->Prev = gCurrNOD->DownstreamPolylines->TailList;
			gCurrNOD->DownstreamPolylines->TailList = gCurrNOD->DownstreamPolylines->TailList->Next;
		}
		// SOMETHING IS WRONG!
		else
		{
			printf("\n Something went wrong in Root NOD assigment at Polyline ID  %d and Compare Polyline ID %d.  No match in x and y to either endpoint. \n STOP. \n", 
				gCurrPLY->ID, ComparePLY->ID);
			return(1);
		}
	}
	
	// Second, now treat the Upstream NOD for the current Polyline...

	// If the HeadList associated with gPEM at cell i_term, j_term is NULL, there are no currently assigned gPEM entities.  
	// We need to create new gPEM entities and create a new upstream NOD here. 
	if((gPEM + i_term * gNumCols + j_term)->HeadList == NULL)
	{
		// Remember: each gPEM cell was allocated as a SEN_VTX_LST with Head and Tail pointers to the VTXs associated 
		// with the gPEM Cell.  Allocate a new HeadList and TaiList for the gPEM - each a pointer to the start and end VTX_LST 
		// associated with this gPEM element
		(gPEM + i_term * gNumCols + j_term)->HeadList = (struct VTX_LST *) calloc(1, sizeof(struct VTX_LST));
		(gPEM + i_term * gNumCols + j_term)->TailList = (struct VTX_LST *) calloc(1, sizeof(struct VTX_LST));	
		
		// Assign the Current PLY's RootVtx to both the HeadList and TailList's VTX
		(gPEM + i_term * gNumCols + j_term)->HeadList->ThisVTX = gCurrPLY->TermVtx;
		
		// Both Head and Tail point to the same VTX intially
		(gPEM + i_term * gNumCols + j_term)->TailList->ThisVTX = gCurrPLY->TermVtx;
		
		// For any first gPEM[i][j] HeadList entry, the previous VTX_LST pointer should be NULL
		// Again, setting the Tail of the List to the same value as the Head initially
		(gPEM + i_term * gNumCols + j_term)->HeadList->Prev = NULL;
		
		(gPEM + i_term * gNumCols + j_term)->TailList->Prev = NULL;
		
		// Set the Next print PEM VTX LST to NULL as well - we will expand the list by expanding the Tail 
		// of the List as needed when this new NOD is later matched as seen in the "else" section below
		(gPEM + i_term * gNumCols + j_term)->HeadList->Next = NULL;

		(gPEM + i_term * gNumCols + j_term)->TailList->Next = NULL;
		
		// Need to assign a new NOD from the allocated space for NODs
		gCurrNOD = (struct NOD *)(gNetworkTable + gNumNODs);
		
		// Assign a VTX to the NOD from current Polyline
		gCurrNOD->ThisVTX= gCurrPLY->TermVtx;
		// Bump the offset index to gNetworkTable for next use
		gNumNODs = gNumNODs + 1;
			
		// Wrt to the current PLY, the current NOD is the Upstream NOD
		// Set Polyline Downstream Node using the next available NOD structure from the Network Table
		gCurrPLY->UpstreamNode = gCurrNOD;
				gCurrNOD = (struct NOD *)(gNetworkTable + gNumNODs);

		// Wrt the current NOD, the current PLY is DOWNSTREAM of the NOD!
		// There may be several PLYs upstream (and/or downstream) of a NOD, therefore we have another 
		// variable-length list, this time for PLYs
		
		// Unlike the gPEM array of cells which were allocated upfront, here we must create the Upstream SEN_PLY_LST for this NOD
		gCurrNOD->DownstreamPolylines = (struct SEN_PLY_LST *)calloc(1, sizeof(struct SEN_PLY_LST));
		
		// Allocate space for the Head and Tail PLY LSTs
		gCurrNOD->DownstreamPolylines->HeadList = (struct PLY_LST *)calloc(1, sizeof(struct PLY_LST));
		gCurrNOD->DownstreamPolylines->TailList = (struct PLY_LST *)calloc(1, sizeof(struct PLY_LST));	
		
		// Set Prev pointers for the Head and Tail LSTs to NULL
		// Initially Head and Tail point to the single entry in the LST
		gCurrNOD->DownstreamPolylines->HeadList->Prev = NULL;
		gCurrNOD->DownstreamPolylines->TailList->Prev = NULL;
		
		// Set Next pointers for the Head and Tail LSTs to NULL
		// Initially Head and Tail point to the single entry in the LST
		gCurrNOD->DownstreamPolylines->HeadList->Next = NULL;
		gCurrNOD->DownstreamPolylines->TailList->Next = NULL;
		
		// Add the current Polyline to the Head and Tail PLY_LSTs for DownstreamPolylines 
		gCurrNOD->DownstreamPolylines->HeadList->ThisPLY = (struct PLY *)gCurrPLY;
		gCurrNOD->DownstreamPolylines->TailList->ThisPLY = (struct PLY *)gCurrPLY;
	}
	// Else, a NOD already exists for this PEM[i][j]. The question becomes, how do we efficiently find this NOD?
	else
	{
		// A NOD already exists. Find the NOD by using a previously assigned entry and it's associated 
		// This VTX and Owner PLY. Use the TailList pointer to quickly find Compare PLY using the ThisVTX Owner PLY. 
		ComparePLY = (gPEM + i_term * gNumCols + j_term)->TailList->ThisVTX->OwnerPLY;

		// Confirm the NOD by comparing x and y of the current Polyline with those of the newly found ComparePLY.
		// An exact check in xy works because coordinates in the SDP data are aligned with the grid post xy's precisely.
		if(gCurrPLY->TermVtx->x == ComparePLY->RootVtx->x && gCurrPLY->TermVtx->y == ComparePLY->RootVtx->y)
		{
			// Polylines match at the ComparePLY Downstream node and the gCurrPLY Upstream Node; for each PLY - 
			// 1) assign PLY Downstream NOD and 2)add current PLY to the list of UPstream PLYs for the NOD
			gCurrNOD = (struct NOD *)ComparePLY->DownstreamNode;
			
			// Wrt to the curent PLY, this is the Upstream NOD
			gCurrPLY->UpstreamNode = gCurrNOD;
			
			// Add the Curr PLY to the NOD's UPSTREAM(!) PLY_LST
			gCurrNOD->ThisVTX = gCurrPLY->RootVtx;
			
			// If we are just adding the second VTX to the PEM cell, special handling is needed for prev and next for the VTX_LSTs
			if((gPEM + i_term * gNumCols + j_term)->TailList->ThisVTX == (gPEM + i_term * gNumCols + j_term)->HeadList->ThisVTX)
			{
				// Initial Tail List update - re-use the existing Tail List
				(gPEM + i_term * gNumCols + j_term)->HeadList->Next = (gPEM + i_term * gNumCols + j_term)->TailList;
				(gPEM + i_term * gNumCols + j_term)->TailList->Prev = (gPEM + i_term * gNumCols + j_term)->HeadList;

				// Assign the current PLY Root VTX to the Tail List This VTX
				(gPEM + i_term * gNumCols + j_term)->TailList->ThisVTX = gCurrPLY->RootVtx;
			}
			// Else we are adding a third, fourth, fifth VTX...always extend from the Tail of the List
			else
			{
				// Update gPEM root TailList (a pointer to a VTX_LST structure)
				(gPEM + i_term * gNumCols + j_term)->TailList->Next = (struct VTX_LST *)calloc(1, sizeof(struct VTX_LST));	
				(gPEM + i_term * gNumCols + j_term)->TailList->Next->Prev = (gPEM + i_term * gNumCols + j_term)->TailList;
				
				// Update to a new List Tail pointer
				(gPEM + i_term * gNumCols + j_term)->TailList = (gPEM + i_term * gNumCols + j_term)->TailList->Next;

				// calloc() NULLs everything implicitly, but be explicit and set Next to NULL 
				(gPEM + i_term * gNumCols + j_term)->TailList->Next = NULL;				

				// Assign the current PLY Root VTX to the Tail List This VTX
				(gPEM + i_term * gNumCols + j_term)->TailList->ThisVTX = gCurrPLY->RootVtx;
			}			
		}
		// Does Current PLY's Root VTX match at the terminal node of Compare PLY? 
		// If so, we have matched one upstream and one downstream PLY at the NOD 
		else if(gCurrPLY->TermVtx->x == ComparePLY->TermVtx->x && gCurrPLY->TermVtx->y == ComparePLY->TermVtx->y)
		{
			// Matches at an upstream node for Compare PLY
			gCurrNOD = (struct NOD *)ComparePLY->UpstreamNode;

			// Note - while the Curr PLY Root VTX matched the Compare PLY Term VTX and it's Upstream NOD,
			// we need to assign the PLY to the NOD's Downstream LST. It may be necessary to draw a picture.
			// Extend the list of polylines and link the LST elements
			gCurrNOD->DownstreamPolylines->TailList->Next = (struct PLY_LST *)calloc(1, sizeof(struct PLY_LST));		
			gCurrNOD->DownstreamPolylines->TailList->Next->Prev = gCurrNOD->DownstreamPolylines->TailList;
			gCurrNOD->DownstreamPolylines->TailList = gCurrNOD->DownstreamPolylines->TailList->Next;
		}
		// SOMETHING IS WRONG!
		else
		{
			printf("\n Something went wrong in Term NOD assigment at Polyline ID  %d and Compare Polyline ID %d.  No match in x and y to either endpoint. \n STOP. \n", 
				gCurrPLY->ID, ComparePLY->ID);
			return(1);
		}
	}
}
 
// Leaving this function in place for now anticipating later inclusion of OWBs and DLDs
// Signed area is a useful test for area orientation (CW, CCW) applicable to SHL of OWBs, DLDs, and CMT
// * Need to amend this function to work off a LST of PLYs or a SHL!
double ComputeClosedLineArea(struct PLY *CurrLine)
{
	double Area = 0.0;
	double TrapezoidArea;		// Partial polygon area for two consecutive points in the boundary
	double yBase = 4408000.0;	// Use the minimum Y for the area of interest as a y-coordinate base to avoid exremely large intermediate trapezoidal areas
	
	struct VTX *FromVtx;
	struct VTX *ToVtx;
	
	FromVtx = (CurrLine)->RootVtx;
	ToVtx = FromVtx->NextVtx;
	
	while(ToVtx != NULL)
	{
		TrapezoidArea = ((ToVtx->y - yBase) + (FromVtx->y - yBase)) * 0.5 * (ToVtx->x - FromVtx->x);
		Area = Area + TrapezoidArea;
		FromVtx = ToVtx;
		ToVtx = ToVtx->NextVtx;
	}
	
	return(Area);
}


int main()
{
	_Bool Status;

	int i, j, idx;
	int Network;
	int NumPLYErrors = 0;
	int NumPEMElements = 0;

	double Length;
	double ElapsedTime = 0.0;
	double CPUTime = 0.0;
	
	char VtxIndex1In[15], VtxIndex2In[15], VtxIndex3In[15];

	struct NOD *Nodes;
	struct LST *ListOfPLYErrors;
	
	clock_t BeginTime, EndTime;

	// Measure execution time
	BeginTime = clock();

	gAllVtxs = (struct VTX *)(calloc(250000, sizeof(struct VTX)));
	
	gAllPolylines = (struct PLY *)(calloc(50000, sizeof(struct PLY)));	
		 
	
	//****************************************************************************
	// I/O Definitions
	// Open input Surface Drainage Paths (SDP) points and DEM as points files, files that
	// create output Network Table and Ordered Surface Drainage Features results
	//****************************************************************************

	gfpSDP_Vtxs = (FILE *)fopen("C://000WORKIndexedTasks 000-999//500 Papers//502 Revisiting DEM Infilling//502_CodeUnitTest//data//RawSDPs.txt","r");
	
	if(gfpSDP_Vtxs == NULL)
	{
		printf("Could not open the Raw SDP Vertices file \n");
		return(1);
	}
	
	gfpDEMAsPoints = (FILE *)fopen("C://000WORKIndexedTasks 000-999//500 Papers//502 Revisiting DEM Infilling//502_CodeUnitTest//data//DEMxyz.txt","r");
	if(gfpDEMAsPoints == NULL)
	{
		printf("Could not open the LiDAR DEM points file. \n");
		return(1);
	}
	gfpSDP_VtxsLogFile = (FILE *)fopen("C://000WORKIndexedTasks 000-999//500 Papers//502 Revisiting DEM Infilling//502_CodeUnitTest//data//LogFile.txt","w");
	
	if(gfpSDP_VtxsLogFile == NULL)
	{
		printf("Could not create the Log file. \n");
		return(1);
	}
	gfpOrderedPlys = (FILE *)fopen("C://000WORKIndexedTasks 000-999//500 Papers//502 Revisiting DEM Infilling//502_CodeUnitTest//data//OrderedPlysFile.txt","w");
	
	if(gfpOrderedPlys == NULL)
	{
		printf("Could not create the output Ordered Surface Drainage PLY Features file. \n");
		return(1);
	}
	
	gfpOrderedVtxs = (FILE *)fopen("C://000WORKIndexedTasks 000-999//500 Papers//502 Revisiting DEM Infilling//502_CodeUnitTest//data//OrderedVtxsFile.txt","w");
	
	if(gfpOrderedVtxs == NULL)
	{
		printf("Could not create the output Ordered Surface Drainage VTX Features file. \n");
		return(1);
	}

	gfpDrainFeaturesNetwork = (FILE *)fopen("C://000WORKIndexedTasks 000-999//500 Papers//502 Revisiting DEM Infilling//502_CodeUnitTest//data//DrainFeasNetwork.txt","w");
	
	if(gfpDrainFeaturesNetwork == NULL)
	{
		printf("Could not create the Drain Lines file. \n");
		return(1);
	}
	fprintf(gfpDrainFeaturesNetwork, "ID Num UpstreamPLYNbrs DownstreamPLYNbrs \n");

	gfpQuededFeatures = (FILE *)fopen("C://000WORKIndexedTasks 000-999//500 Papers//502 Revisiting DEM Infilling//502_CodeUnitTest//data//QueuedEdit.txt","w");
	
	if(gfpQuededFeatures == NULL)
	{
		printf("Could not create the output Queued Edit Features file. \n");
		return(1);
	}
	
	// Initialize the memory pool for 250K gAllVtxs (202145 0bserved in test data)
	// This limit may need to be raised when densification of the SDP vertices occurs (500K)
	for(idx = 0; idx < 250000; idx++)
	{
		(gAllVtxs+idx)->x = 0.0;
		(gAllVtxs+idx)->y = 0.0;
		(gAllVtxs+idx)->z = 0.0;			
		(gAllVtxs+idx)->NormDist = 0.0;
		(gAllVtxs+idx)->OwnerPLY = NULL;
		(gAllVtxs+idx)->AccumUpslopeArea = -1.0;
		(gAllVtxs+idx)->NextVtx = NULL;
		(gAllVtxs+idx)->PrevVtx = NULL;
	}
	
	// Initialize the memory pool for SDP gAllPolylines (38330 observed in test data)
	for(idx=0; idx < 50000; idx++)
	{
		(gAllPolylines + idx)->ID = idx;
		(gAllPolylines + idx)->RootVtx = NULL;
		(gAllPolylines + idx)->TermVtx = NULL;
		(gAllPolylines + idx)->DownstreamNode = NULL;
		(gAllPolylines + idx)->UpstreamNode = NULL;
		(gAllPolylines + idx)->Length = 0.0;
		(gAllPolylines + idx)->OwnerSHL = NULL;
		(gAllPolylines + idx)->EntityLeft = NULL;
		(gAllPolylines + idx)->EntityRight = NULL;
	}
	
	//****************************************************************************
	// Load the DEM array
	//****************************************************************************
	gNumPosts = LoadDEMPosts();
	
	//***************************************************************************
	// Load the Surface Drainage Path (SDP) vertices
	// Sets values for gNumVTXs, gNumPLYs
	//***************************************************************************
	LoadSDPVtxs();

	// Report results of vertices and polyline read to output
	printf("\n %d Vertices read and processed for %d polylines. \n", gNumVTXs, gNumPLYs);
	
	//Gather statistics on difference in elevation between Polyline Root and Term VTXs
	//Evaluate gross monotonicity and Threshold Diff used later as a sanity check for whether or not to reverse "flat" or "inverted streams" 
	EvaluateThresholdElevDiff();
	
	printf(" \n Endpoint elevation differences successfully analyzed and a Threshold value of %lf was developed. \n", gThresholdDiff);
	
	//*******************************************************************************************************
	// Set PEM pointer entity type to POE and initialize the PEM LST entity and prev and next pointers to NULL
	// The entity for this LST is a POE and the PEM[i][j]->Entity is NULLed -  we must allocate the POE as we go
	//*******************************************************************************************************
	NumPEMElements = CreatePEM();
	
	//*******************************************************************************************************
	// Now create the Nodes to Polylines LST table. Polyline LSTs are already established,we will just need to 
	// extend the LST and set the appropriate values.
	//*******************************************************************************************************
	CreateNODTable();
	
	// Note switched from jdx to idx as the index to Polylines
	idx = 0;
	
	// Now map the Polyline root and term VTX pointers to a PEM grid cell
	while(idx < gNumPLYs)
	{
		// Workhorse routine - most of the work of the program is done here in a single pass through the Polylines!
		// Basic polyline to node checks, including a monotonic check based on polyline endpoint Zs and a second pass through
		// all NODs comparing opposing endpoints of polyline features matched ata common NOD
		// When complete:
		//  * all PLYs and NODs are connected  
		//  * the PEM elements containing drain vertices have at least one POE entry (more for PEM elements at NOD locations)
		//  * the Surface Drainage Polyline (SDP) Surface Drainage Features Network Table (SDN) is complete, revealing the lists of upstream and 
		//    downstream polylines at the nodes
		//  * the initial version of the SDN is complete and the geometry is ordered and ready for densification if needed
		MapNODs(idx);
		
		// Flush the output buffers to ensure availability of interim data
		fflush(gfpOrderedPlys);
		fflush(gfpOrderedVtxs);
		
		OutputWorkingtNetworkTable(idx);
		fflush(gfpDrainFeaturesNetwork);
		
		idx = idx + 1;
	}
	


		
	// One last check for global monotonicity of the Network now
	NumPLYErrors = PiecewiseValidationOfGlobalNetworkConnectivity();
	
	if(NumPLYErrors > 0)
	{
		printf("\n Overall monotonicity of the network failed with %OrderedPLd errors.  See the output file QueuedEdit.txt for details. \n", NumPLYErrors);
	}
	else
	{
		printf("\n Congratulations! Overall monotonicity of the network has been validated by a cumulative piecewise procedure. \n");
	}
		
	// Output the Edit Queue and the ambiguous drains results
	UpdateEditQueueForNetwork();

	printf (" \n %d treated SDP polylines created. \n", gNumPLYs);
 	printf (" \n %d treated polylines with monotonicity errors. \n", gNumPLYErrors); 
	
	Status = OutputAll();
	
	if(Status == FALSE)
	{
		printf("\n Failed to output all data. \n STOP.\n");
		return(1);
	}
	else
	{
		printf(" \n Program 502ProcessB_TopoMono_Checks.c successfully completed. \n");
		
		EndTime = clock();

		CPUTime = (double)(EndTime - BeginTime);
		ElapsedTime = (double) CPUTime / CLOCKS_PER_SEC;
		
		printf("\n Executed in %lf seconds wall time and %lf CPU seconds", ElapsedTime, CPUTime);
	}
};
