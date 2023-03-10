#ifndef FILE_VISPAR
#define FILE_VISPAR

class VisualizationParameters
{
public:
  double lightamb;
  double lightdiff;
  double lightspec;
  double shininess;
  double transp;
  int locviewer;
  char selectvisual[20];
  int showstltrias;
  
  Vec3d clipnormal;
  double clipdist;
  int clipenable;
  int clipplanetimestamp;

  int colormeshsize;

  int drawfilledtrigs;
  int drawbadels;
  int drawoutline;
  int drawedges;
  int subdivisions;

  int drawprisms;
  int drawpyramids;
  int drawhexes;
  double shrink;
  int drawidentified;
  int drawpointnumbers;
  int drawedgenumbers;
  int drawfacenumbers;
  int drawelementnumbers;
  int drawdomainsurf;
  int drawtets;
  int drawtetsdomain;

  int drawededges;
  int drawedpoints;
  int drawedpointnrs;
  int drawedtangents;
  int drawededgenrs;
  int drawmetispartition;

  int drawcurveproj;
  int drawcurveprojedge;
  

  int centerpoint;
  int drawelement;

  // stl:
  int stlshowtrias;
  int stlshowfilledtrias;
  int stlshowedges;
  int stlshowmarktrias;
  int stlshowactivechart;
  int stlchartnumber;
  int stlchartnumberoffset;

  // occ:
  int occshowvolumenr;
  bool occshowsurfaces;
  bool occshowedges;
  bool occvisproblemfaces;
  bool occzoomtohighlightedentity;
  double occdeflection;

  bool whitebackground;
  int stereo;
  bool usedispllists;
  bool drawcoordinatecross;
  bool drawcolorbar;
  bool drawnetgenlogo;

  bool use_center_coords;
  double centerx,centery,centerz;

  
public:
  VisualizationParameters();
};
extern VisualizationParameters vispar;

#endif
