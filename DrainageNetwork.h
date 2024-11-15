#ifndef HydroFeatures
#define HydroFeatures 1

struct CMT									// A catchment feature represented alternatively as a raster and collection of grid cells or a vector as a closed Boundary
{
	unsigned long Prop;						// Property bits identifying the Type of catchment and processing flags
	int ID;
	int MinRow;								// Raster representation of the CMT extent
	int MinCol;
	int MaxRow;
	int MaxCol;
	struct PLY_LST *CatchBDY;				// A collection of all Polylines on the CMT Boundary (vector-topological representation)	
	struct VTX_LST *PourPoint;				// A collection of all Pour Points (NODs) on the Catchment Boundary
};

struct DEP									// A depression feature
{
	int MinRow;								// MBR Raster representation of the DEP
	int MinCol;
	int MaxRow;
	int MaxCol;
	int Prop;								// Property bits identifying the Type of catchment and processing flags
	int ID;
	struct NOD *MinPoint;					// The lowest point in the depressed area, centermost if the minimum point occurs in a flat area
	struct OWB *DEPAsOWB;					// A pointer to an OWB feature if one has been previously established
};

struct DLD									// A double-line drain feature
{
	int Prop;								// Property bits for processing and state of the DLD (raster or vector)
	int ID;									// A sequential ID number for this DLD
	struct CMT *WatershedBDY;				// An ordered collection of all Polylines on the DLD Boundary extent(vector-topological representation)	
	struct PLY_LST *RightSHL;				// An ordered collection of Polylines making up the right shoreline component of the DLD wrt downstream direction
	struct PLY_LST *LeftSHL;				// A collection of Polylines making up the left shoreline component of the DLD wrt downstream direction
	struct CMT *CatchmentBdy;				// Allows the alternative representations of a DLD
};

struct HRU									// An atomic hydrologic response unit of analysis in hydrologic modelling homogenous with respect to hillslope position, slope land cover and soils
{											// A subdivision of a path of overland flow
	struct POF *ParentPOF;					// The parent POF for this HRU
	int Cover;								// Land Cover index to hydrologic characteristics
	int Hillslope;							// Hillslope position (e.g., Summital Convexity, Free face, Rectilinear, Basal Concavity - see https://ebooks.inflibnet.ac.in/geop11/chapter/slope-elements-and-slope-evolution/
	double Slope;							// Prevailing slope over HRU
	int Soil;								// Soils index to hydrologic characteristics
	double DepthWater;						// Estimated depth to subsurface water flow 
};

struct ISL									// An island feature
{
	int Prop;								// Property bits for processing and state of the ISL (raster or vector)
	int Type;								// A type code that may be needed to distinguish Island types
	struct SHL *IslandSHL;					// A list of shoreline PLYs that make up the CCW island outer boundary
};

struct NOD									// A topological node connecting Polylines and representing important point features
{
	int ID;
	struct VTX *ThisVTX;						// 
	struct SEN_PLY_LST *DownstreamPolylines;	// A collection of pointers to Polylines immediately downstream of the node
	struct SEN_PLY_LST *UpstreamPolylines;		// A collection of pointers to Polylines immediately upstream of the node
};

struct SEN_NOD_LST							// Holds sentinnel pointers to the head and tail of the NOD_LST
{
	struct PLY_LST *HeadList;
	struct PLY_LST *TailList;
};

struct OWB									// An open water body feature
{
	int Prop;								// Processsing flags for the OWB
	struct SHL *Boundary;					// The shoreline Boundary including islands
	struct CMT *CatchmentBdy;				// A pointer to the grid cell collection
	struct VTX_LST *PourPoints;				// A collection of pointers to vertices that are "saddle points" and potential pour points on the catchment Boundary
};

struct PLY
{
	// NOTE - don't you forget it: ideally, root = downstream and term = upstream.  The vertices in lines always proceed from root to term.  
	// By convention, that progression should also be low to high Z. If TermVtx->z - RootVtx->z is negative, the Polyline is 
	// "flowing uphill" and should be reversed.
	int ID;
	int NetworkIndex;
	double Length;
	struct VTX *RootVtx;
	struct VTX *TermVtx;
	struct NOD *DownstreamNode;
	struct NOD *UpstreamNode;
	struct SHL *OwnerSHL;
	struct POF *EntityLeft;
	struct POF *EntityRight;
};

struct PLY_LST
{
	struct PLY *ThisPLY;
	struct PLY_LST *Next;
	struct PLY_LST *Prev;
};

struct SEN_PLY_LST							// Holds sentinnel pointers to the head and tail of the PLY_LST
{
	struct PLY_LST *HeadList;
	struct PLY_LST *TailList;
};


struct POE
{
	struct VTX_LST *PEMVtx;					// A pointer to the first of possibly multiple vertices (ie, at a NOD) in the PEM
};

struct POF
{
	int Prop;								// Property bits for processing and state of the DLD (raster or vector)
	int ID;
	struct SHL *WaterBDY;					// Or a collection of all Polylines on the DLD Boundary extent(vector-topological representation)	
	struct PLY_LST *POF_RHC;				// A collection of Polylines making up the right hand component of the POF
	struct PLY_LST *POF_LHC;				// A collection of Polylines making up the left hand component of the POF
	struct CMT *CatchmentBdy;				// Allows the alternative representations of a DLD
};

struct PPT
{
	struct VTX *PourPoint;	
};

struct RCH
{
	struct POF *ParentPath;
	struct PLY *Feature1;
	double Dist1;
	struct PLY *Feature2;
	double Dist2;
};

struct SDP									// A surface drainage path - a centerline for the channel
{
	int Prop;								// Processsing flags for the SDP
	struct LST *Centerline;					// A collection of Polylines making up the centerline representation of the SDP
	struct CMT *CatchBdy;
	struct DLD *SDPAsDLD;					// A pointer to a DLD feature if the SDP has been promoted to a DLD
};

struct SHL									// Shoreline feature component - the precise location of the land and water interface
{
	int Prop;								// Processsing flags for the SHL
	void *OwnerFea;							// An undifferentiated pointer to a feature of this shoreline - a SDP centerline, a OWB, a DLD, or a POF(?) 
	struct PLY_LST *OuterBdy;				// A collection of all Polylines making up the shoreline outer extent: orientation is CW 
	struct PLY_LST *Islands;				// A collection of simple CCW-closed boundaries making up islands interior to the shoreline
};

struct VTX
{
	double x;
	double y;
	double z;
	double NormDist;
	double AccumUpslopeArea;
	struct PLY *OwnerPLY;
	struct VTX *NextVtx;
	struct VTX *PrevVtx;
};

struct VTX_LST
{
	struct VTX *ThisVTX;
	struct VTX_LST *Next;
	struct VTX_LST *Prev;
};

struct SEN_VTX_LST							// Holds sentinnel pointers to the head and tail of the VTX_LST
{
	struct VTX_LST *HeadList;
	struct VTX_LST *TailList;
};
#endif