(*^

::[	Information =

	"This is a Mathematica Notebook file.  It contains ASCII text, and can be
	transferred by email, ftp, or other text-file transfer utility.  It should
	be read or edited using a copy of Mathematica or MathReader.  If you 
	received this as email, use your mail application or copy/paste to save 
	everything from the line containing (*^ down to the line containing ^*)
	into a plain text file.  On some systems you may have to give the file a 
	name ending with ".ma" to allow Mathematica to recognize it as a Notebook.
	The line below identifies what version of Mathematica created this file,
	but it can be opened using any other version as well.";

	FrontEndVersion = "NeXT Mathematica Notebook Front End Version 2.2";

	NeXTStandardFontEncoding; 
	
	fontset = title, inactive, noPageBreakBelow, noPageBreakInGroup, nohscroll, preserveAspect, groupLikeTitle, center, M7, bold, L1, e8,  24, "Times"; ;
	fontset = subtitle, inactive, noPageBreakBelow, noPageBreakInGroup, nohscroll, preserveAspect, groupLikeTitle, center, M7, bold, L1, e6,  18, "Times"; ;
	fontset = subsubtitle, inactive, noPageBreakBelow, noPageBreakInGroup, nohscroll, preserveAspect, groupLikeTitle, center, M7, italic, L1, e6,  14, "Times"; ;
	fontset = section, inactive, noPageBreakBelow, nohscroll, preserveAspect, groupLikeSection, grayBox, M22, bold, L1, a20,  18, "Times"; ;
	fontset = subsection, inactive, noPageBreakBelow, nohscroll, preserveAspect, groupLikeSection, blackBox, M19, bold, L1, a15,  14, "Times"; ;
	fontset = subsubsection, inactive, noPageBreakBelow, nohscroll, preserveAspect, groupLikeSection, whiteBox, M18, bold, L1, a12,  12, "Times"; ;
	fontset = text, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7, L1,  12;
	fontset = smalltext, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7, L1,  10, "Times"; ;
	fontset = input, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeInput, M42, N23, bold, L1,  12, "Courier"; ;
	fontset = output, output, inactive, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeOutput, M42, N23, L-5,  12, "Courier"; ;
	fontset = message, inactive, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeOutput, M42, N23, L1,  12, "Courier"; ;
	fontset = print, inactive, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeOutput, M42, N23, L1,  12, "Courier"; ;
	fontset = info, inactive, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeOutput, M42, N23, L1,  12, "Courier"; ;
	fontset = postscript, PostScript, formatAsPostScript, output, inactive, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeGraphics, M7, l34, w282, h287, L1,  12, "Courier"; ;
	fontset = name, inactive, noPageBreakInGroup, nohscroll, preserveAspect, M7, italic, B65535, L1,  10, "Times"; ;
	fontset = header, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7, italic, L1,  12, "Times"; ;
	fontset = leftheader,  12;
	fontset = footer, inactive, nohscroll, noKeepOnOnePage, preserveAspect, center, M7, italic, L1,  12, "Times"; ;
	fontset = leftfooter,  12;
	fontset = help, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7, L1,  12;
	fontset = clipboard, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7, L1,  12;
	fontset = completions, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7, L1,  12, "Courier"; ;
	fontset = special1, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7, L1,  12;
	fontset = special2, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7, L1,  12;
	fontset = special3, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7, L1,  12;
	fontset = special4, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7, L1,  12;
	fontset = special5, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7, L1,  12;
	currentKernel; 
]
:[font = title; inactive; preserveAspect]
PolytopeSkeleton.m
Visualizing Convex Polytope Skeletons (with cddmathlink)
:[font = subtitle; inactive; preserveAspect]
Komei Fukuda, fukuda@ifor.math.ethz.ch
Swiss Federal Institute of Technology, Lausanne and Zurich
January 9, 2000
;[s]
1:0,0;113,-1;
1:1,15,12,Times,3,17,0,0,0;
:[font = section; inactive; initialization; Cclosed; preserveAspect; startGroup]
Preparation (reading Packages)
:[font = input; initialization; preserveAspect]
*)
Off[General::spell1]; Off[General::spell]
(*
:[font = input; initialization; preserveAspect]
*)
$Path = Append[$Path,"~/Math"]; 
(*
:[font = text; inactive; initialization; preserveAspect]
We use an extra graphics package available from MathSource.    It is called View3D which is a part of ExtendGraphics  package wrtten by  Tom Wickham-Jones .  The package is available at http://www.mathsource.com/.  If you install it in a directory which in not in the $Path, you need to append the directory to $Path.
;[s]
5:0,0;48,1;58,2;186,3;212,4;318,-1;
5:1,11,8,Times,0,12,0,0,0;1,10,8,Times,2,12,0,0,0;1,11,8,Times,0,12,0,0,0;1,10,8,Courier,1,12,0,0,0;1,11,8,Times,0,12,0,0,0;
:[font = input; initialization; preserveAspect; startGroup]
*)
$Path
(*
:[font = output; output; inactive; initialization; preserveAspect; endGroup]
{".", "~", "~/Library/Mathematica/Packages", "/LocalApps/Mathematica22.i3\
 
   86.app/Library/Mathematica/Packages", "/LocalLibrary/Mathematica/Packa\
 
   ges", "/LocalApps/Mathematica22.i386.app/Install/Preload", "/LocalApps\
 
   /Mathematica22.i386.app/StartUp", "~/Math", "~/Math", "~/Math"}
;[o]
{., ~, ~/Library/Mathematica/Packages, 
 
  /LocalApps/Mathematica22.i386.app/Library/Mathematica/Packages, 
 
  /LocalLibrary/Mathematica/Packages, 
 
  /LocalApps/Mathematica22.i386.app/Install/Preload, 
 
  /LocalApps/Mathematica22.i386.app/StartUp, ~/Math, ~/Math, ~/Math}
:[font = input; initialization; preserveAspect]
*)
Needs["ExtendGraphics`View3D`"];
(*
:[font = text; inactive; initialization; preserveAspect]
cddmathlink is a MathLink  version of libcdd for doing the vertex enumeration and the facet enumeration of a convex polyhedron.  cddlib-080.tar.gz is available from Fukuda's Homepage:   http:///www.ifor.math.ethz.ch/staff/fukuda/fukuda.html  .   Look for "cdd Homepage".
;[s]
3:0,0;17,1;25,2;270,-1;
3:1,11,8,Times,0,12,0,0,0;1,10,8,Times,2,12,0,0,0;1,11,8,Times,0,12,0,0,0;
:[font = input; initialization; preserveAspect; startGroup]
*)
cddml=Install["~/Math/cddmathlink"]
(*
:[font = output; output; inactive; initialization; preserveAspect; endGroup]
LinkObject["~/Math/cddmathlink", 22, 22]
;[o]
LinkObject[~/Math/cddmathlink, 22, 22]
:[font = text; inactive; initialization; preserveAspect]
The following packages come with libcdd package and should be installed in your favorite directory.  Here I assume it is in "~/Math" directory.   They are separately available from Fukuda's Homepage: http:///www.ifor.math.ethz.ch/staff/fukuda/fukuda.html  .  Look for Mathematica Projects.
;[s]
3:0,0;268,1;279,2;289,-1;
3:1,11,8,Times,0,12,0,0,0;1,10,8,Times,2,12,0,0,0;1,11,8,Times,0,12,0,0,0;
:[font = input; initialization; preserveAspect]
*)
<< PolytopeSkeleton.m
(*
:[font = input; initialization; preserveAspect]
*)
<< IOPolyhedra.m
(*
:[font = input; initialization; preserveAspect; endGroup]
*)
<< UnfoldPolytope.m
(*
:[font = section; inactive; Cclosed; preserveAspect; startGroup]
Compute theVertices of a Polytope by using cddmathlink
:[font = input; preserveAspect; startGroup]
FileNames["*","~/Math/ine_3d"]
:[font = output; output; inactive; preserveAspect; endGroup]
{"~/Math/ine_3d/cube3.ine", "~/Math/ine_3d/cubocta.ine", "~/Math/ine_3d/c\
 
   uboctaT.ine", "~/Math/ine_3d/dodeca.ine", "~/Math/ine_3d/dodecaT.ine"\
 
   , "~/Math/ine_3d/grcubocta.ine", "~/Math/ine_3d/grcuboctaT.ine", "~/Ma\
 
   th/ine_3d/hexocta.ine", "~/Math/ine_3d/hexoctaT.ine", "~/Math/ine_3d/i\
 
   cododeca.ine", "~/Math/ine_3d/icododecaT.ine", "~/Math/ine_3d/rcubocta\
 
   .ine", "~/Math/ine_3d/rcuboctaT.ine", "~/Math/ine_3d/rhomtria.ine"\
 
   , "~/Math/ine_3d/rhomtriaT.ine"}
;[o]
{~/Math/ine_3d/cube3.ine, ~/Math/ine_3d/cubocta.ine, 
 
  ~/Math/ine_3d/cuboctaT.ine, ~/Math/ine_3d/dodeca.ine, 
 
  ~/Math/ine_3d/dodecaT.ine, ~/Math/ine_3d/grcubocta.ine, 
 
  ~/Math/ine_3d/grcuboctaT.ine, ~/Math/ine_3d/hexocta.ine, 
 
  ~/Math/ine_3d/hexoctaT.ine, ~/Math/ine_3d/icododeca.ine, 
 
  ~/Math/ine_3d/icododecaT.ine, ~/Math/ine_3d/rcubocta.ine, 
 
  ~/Math/ine_3d/rcuboctaT.ine, ~/Math/ine_3d/rhomtria.ine, 
 
  ~/Math/ine_3d/rhomtriaT.ine}
:[font = text; inactive; preserveAspect]
The file test.ine is the main output file of cdd+, which gives the inequality (facet) representation of the polytope.
:[font = input; preserveAspect; startGroup]
MatrixForm[inedata=
	ReadPolyhedraData["~/Math/ine_3d/hexocta.ine"]];
:[font = print; inactive; preserveAspect; endGroup]
m=48,   n=4 type=integer
:[font = input; preserveAspect]
amat=-Transpose[Drop[Transpose[inedata],1]];
bvec=Transpose[inedata][[1]];
:[font = input; preserveAspect; startGroup]
{m,d}=Dimensions[inedata]
:[font = output; output; inactive; preserveAspect; endGroup]
{48, 4}
;[o]
{48, 4}
:[font = input; preserveAspect]
{{extlist,linearity},ecdlist,eadlist}=
	AllVerticesWithAdjacency[m,d,Flatten[inedata]];
:[font = input; preserveAspect; startGroup]
vlist = Map[Drop[#, 1]&, extlist]
:[font = output; output; inactive; preserveAspect; endGroup]
{{-0.60000000000000009, 0, 0.60000000000000018}, {0, 0, 1.}, 
 
  {-0.50000000000000018, 0.49999999999999956, 0.50000000000000009}, 
 
  {0, 0, -1.}, {-0.50000000000000009, -0.5, -0.49999999999999991}, 
 
  {0.60000000000000044, -0.59999999999999982, 0}, 
 
  {0.50000000000000044, 0.49999999999999964, -0.49999999999999991}, 
 
  {0, 0.59999999999999973, 0.60000000000000062}, 
 
  {-0.60000000000000089, 0.59999999999999911, 0}, 
 
  {0, -0.60000000000000009, 0.60000000000000009}, 
 
  {-0.59999999999999982, 0, -0.59999999999999982}, 
 
  {0.60000000000000062, 0.59999999999999956, 0}, 
 
  {0.60000000000000053, 0, -0.6}, {0, 1., 0}, 
 
  {-0.49999999999999991, 0.5000000000000008, -0.49999999999999964}, 
 
  {0, 0.59999999999999982, -0.60000000000000009}, 
 
  {0.49999999999999982, -0.50000000000000018, 0.50000000000000018}, 
 
  {0.49999999999999991, 0.49999999999999991, 0.50000000000000053}, 
 
  {0.59999999999999973, 0, 0.60000000000000036}, {1., 0, 0}, 
 
  {-0.60000000000000018, -0.59999999999999964, 0}, {-1., 0, 0}, 
 
  {-0.50000000000000018, -0.50000000000000009, 0.49999999999999991}, 
 
  {0.50000000000000071, -0.49999999999999858, -0.50000000000000009}, 
 
  {0, -1., 0}, {0, -0.59999999999999982, -0.60000000000000009}}
;[o]
{{-0.60000000000000009, 0, 0.60000000000000018}, {0, 0, 1.}, 
 
  {-0.50000000000000018, 0.49999999999999956, 0.50000000000000009}, 
 
  {0, 0, -1.}, {-0.50000000000000009, -0.5, -0.49999999999999991}, 
 
  {0.60000000000000044, -0.59999999999999982, 0}, 
 
  {0.50000000000000044, 0.49999999999999964, -0.49999999999999991}, 
 
  {0, 0.59999999999999973, 0.60000000000000062}, 
 
  {-0.60000000000000089, 0.59999999999999911, 0}, 
 
  {0, -0.60000000000000009, 0.60000000000000009}, 
 
  {-0.59999999999999982, 0, -0.59999999999999982}, 
 
  {0.60000000000000062, 0.59999999999999956, 0}, 
 
  {0.60000000000000053, 0, -0.6}, {0, 1., 0}, 
 
  {-0.49999999999999991, 0.5000000000000008, -0.49999999999999964}, 
 
  {0, 0.59999999999999982, -0.60000000000000009}, 
 
  {0.49999999999999982, -0.50000000000000018, 0.50000000000000018}, 
 
  {0.49999999999999991, 0.49999999999999991, 0.50000000000000053}, 
 
  {0.59999999999999973, 0, 0.60000000000000036}, {1., 0, 0}, 
 
  {-0.60000000000000018, -0.59999999999999964, 0}, {-1., 0, 0}, 
 
  {-0.50000000000000018, -0.50000000000000009, 0.49999999999999991}, 
 
  {0.50000000000000071, -0.49999999999999858, -0.50000000000000009}, 
 
  {0, -1., 0}, {0, -0.59999999999999982, -0.60000000000000009}}
:[font = input; preserveAspect]
invertpos[l_List, j_]:=
	Module[{i,pos={}},
		Do[
			If[Position[l[[i]],j]!={},AppendTo[pos,i]],
			{i,1,Length[l]}
		];
		pos
	];
dualizeIndices[ecd_List, ineLen_Integer]:=
	Map[invertpos[ecd,#]&, Range[ineLen]];
:[font = input; preserveAspect; startGroup]
icdlist=dualizeIndices[ecdlist, Length[inedata]]
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
{{17, 19, 20}, {2, 18, 19}, {2, 17, 19}, {18, 19, 20}, {6, 17, 25}, 
 
  {6, 24, 25}, {6, 17, 20}, {6, 20, 24}, {14, 15, 16}, {4, 15, 16}, 
 
  {7, 14, 16}, {4, 7, 16}, {4, 5, 26}, {5, 25, 26}, {4, 24, 26}, 
 
  {24, 25, 26}, {1, 22, 23}, {1, 2, 3}, {1, 3, 22}, {1, 2, 23}, 
 
  {21, 23, 25}, {21, 22, 23}, {5, 21, 25}, {5, 21, 22}, {4, 5, 11}, 
 
  {11, 15, 22}, {4, 11, 15}, {5, 11, 22}, {9, 14, 15}, {3, 9, 14}, 
 
  {9, 15, 22}, {3, 9, 22}, {10, 17, 25}, {2, 10, 17}, {10, 23, 25}, 
 
  {2, 10, 23}, {2, 8, 18}, {8, 14, 18}, {2, 3, 8}, {3, 8, 14}, 
 
  {7, 13, 20}, {4, 13, 24}, {4, 7, 13}, {13, 20, 24}, {7, 12, 14}, 
 
  {12, 14, 18}, {7, 12, 20}, {12, 18, 20}}
;[o]
{{17, 19, 20}, {2, 18, 19}, {2, 17, 19}, {18, 19, 20}, {6, 17, 25}, 
 
  {6, 24, 25}, {6, 17, 20}, {6, 20, 24}, {14, 15, 16}, {4, 15, 16}, 
 
  {7, 14, 16}, {4, 7, 16}, {4, 5, 26}, {5, 25, 26}, {4, 24, 26}, 
 
  {24, 25, 26}, {1, 22, 23}, {1, 2, 3}, {1, 3, 22}, {1, 2, 23}, 
 
  {21, 23, 25}, {21, 22, 23}, {5, 21, 25}, {5, 21, 22}, {4, 5, 11}, 
 
  {11, 15, 22}, {4, 11, 15}, {5, 11, 22}, {9, 14, 15}, {3, 9, 14}, 
 
  {9, 15, 22}, {3, 9, 22}, {10, 17, 25}, {2, 10, 17}, {10, 23, 25}, 
 
  {2, 10, 23}, {2, 8, 18}, {8, 14, 18}, {2, 3, 8}, {3, 8, 14}, 
 
  {7, 13, 20}, {4, 13, 24}, {4, 7, 13}, {13, 20, 24}, {7, 12, 14}, 
 
  {12, 14, 18}, {7, 12, 20}, {12, 18, 20}}
:[font = section; inactive; Cclosed; preserveAspect; startGroup]
Drawing the Skeleton of a Polytope
:[font = subsection; inactive; Cclosed; preserveAspect; startGroup]
Defining 3D objects.  
:[font = input; preserveAspect; startGroup]
skel3D[vp_]:= 
	Graphics3D[
		VisibleSkeleton[vlist, ecdlist, eadlist, 
     	{amat, bvec}, vp]
     ]; 

:[font = input; preserveAspect; endGroup; endGroup]
facets=Map[(Part[vlist,#]) &, icdlist];
facets1=OrderFacets[facets];
solid3D:= Graphics3D[Polygon /@ facets1]; 

:[font = subsection; inactive; Cclosed; preserveAspect; startGroup]
Draw a Polytope  by Built-In Graphics
:[font = input; preserveAspect; startGroup; animationSpeed = 29]
Do[ uvp=circle[a,5,2];
	Show[solid3D, 
		Boxed->False, 
 		ViewPoint -> getMmaViewPoint[uvp,skel3D[uvp]],
 		SphericalRegion->True
 	],
 	{a,Pi/11,Pi,Pi/18}
 ]
:[font = postscript; PostScript; formatAsPostScript; output; inactive; preserveAspect; pictureLeft = 34; pictureWidth = 282; pictureHeight = 282; endGroup; endGroup]
%!
%%Creator: Mathematica
%%AspectRatio: 1 
MathPictureStart
%% Graphics3D
/Courier findfont 10  scalefont  setfont
% Scaling calculations
0.02381 0.952381 0.02381 0.952381 [
[ 0 0 0 0 ]
[ 1 1 0 0 ]
] MathScale
% Start of Graphics
1 setlinecap
1 setlinejoin
newpath
[ ] 0 setdash
0 g
0 0 m
1 0 L
1 1 L
0 1 L
closepath
clip
newpath
p
P
p
.002 w
.64938 .64203 m .67152 .53478 L .55958 .57536 L closepath p
.334 .121 .435 r
F P
s
P
p
.002 w
.63995 .42369 m .67152 .53478 L .55958 .57536 L closepath p
.434 .359 .673 r
F P
s
P
p
.002 w
.42683 .44953 m .53652 .42583 L .55958 .57536 L closepath p
.713 .619 .75 r
F P
s
P
p
.002 w
.53652 .42583 m .55958 .57536 L .63995 .42369 L closepath p
.561 .511 .757 r
F P
s
P
p
.002 w
.41118 .56041 m .42683 .44953 L .55958 .57536 L closepath p
.774 .613 .68 r
F P
s
P
p
.002 w
.42213 .65964 m .41118 .56041 L .55958 .57536 L closepath p
.765 .483 .507 r
F P
s
P
p
.002 w
.53938 .67959 m .64938 .64203 L .55958 .57536 L closepath p
.417 .094 .308 r
F P
s
P
p
.002 w
.53938 .67959 m .42213 .65964 L .55958 .57536 L closepath p
.66 .311 .37 r
F P
s
P
p
.002 w
.5 .28594 m .53652 .42583 L .42683 .44953 L closepath p
.721 .705 .835 r
F P
s
P
p
.002 w
.27448 .52459 m .41118 .56041 L .42683 .44953 L closepath p
.861 .697 .677 r
F P
s
P
p
.002 w
.36746 .3862 m .27448 .52459 L .42683 .44953 L closepath p
.888 .809 .765 r
F P
s
P
p
.002 w
.36746 .3862 m .5 .28594 L .42683 .44953 L closepath p
.806 .811 .855 r
F P
s
P
p
.002 w
.5 .28594 m .53652 .42583 L .63995 .42369 L closepath p
.547 .593 .851 r
F P
s
P
p
.002 w
.41118 .56041 m .27448 .52459 L .42213 .65964 L closepath p
.863 .558 .47 r
F P
s
P
p
.002 w
.74854 .4729 m .67152 .53478 L .63995 .42369 L closepath p
.175 .194 .64 r
F P
s
P
p
.002 w
.64017 .34908 m .5 .28594 L .63995 .42369 L closepath p
.342 .523 .889 r
F P
s
P
p
.002 w
.64017 .34908 m .74854 .4729 L .63995 .42369 L closepath p
.061 .255 .748 r
F P
s
P
p
.002 w
.5 .74578 m .35636 .65466 L .42213 .65964 L closepath p
.626 .065 0 r
F P
s
P
p
.002 w
.27448 .52459 m .35636 .65466 L .42213 .65964 L closepath p
.862 .427 .174 r
F P
s
P
p
.002 w
.53938 .67959 m .5 .74578 L .42213 .65964 L closepath p
.566 .066 .028 r
F P
s
P
p
.002 w
.27448 .52459 m .36746 .3862 L .3433 .351 L closepath p
.918 .978 .732 r
F P
s
P
p
.002 w
.36746 .3862 m .3433 .351 L .5 .28594 L closepath p
.79 .986 .894 r
F P
s
P
p
.002 w
.67152 .53478 m .64938 .64203 L .74854 .4729 L closepath p
0 0 .253 r
F P
s
P
p
.002 w
.5 .28594 m .45578 .29834 L .3433 .351 L closepath p
.361 .819 .715 r
F P
s
P
p
.002 w
.64017 .34908 m .59005 .31539 L .5 .28594 L closepath p
0 .362 .749 r
F P
s
P
p
.002 w
.59005 .31539 m .45578 .29834 L .5 .28594 L closepath p
.255 0 0 r
F P
s
P
p
.002 w
.53938 .67959 m .5 .74578 L .64938 .64203 L closepath p
.15 0 0 r
F P
s
P
p
.002 w
.33138 .59194 m .35636 .65466 L .27448 .52459 L closepath p
0 0 .487 r
F P
s
P
p
.002 w
.30262 .45997 m .3433 .351 L .27448 .52459 L closepath p
0 0 0 r
F P
s
P
p
.002 w
.30262 .45997 m .33138 .59194 L .27448 .52459 L closepath p
0 0 .455 r
F P
s
P
p
.002 w
.74854 .4729 m .65265 .63107 L .64938 .64203 L closepath p
.644 .823 .403 r
F P
s
P
p
.002 w
.5 .74578 m .65265 .63107 L .64938 .64203 L closepath p
.514 .877 .654 r
F P
s
P
p
.002 w
.64017 .34908 m .59005 .31539 L .74854 .4729 L closepath p
.671 .159 0 r
F P
s
P
p
.002 w
.33138 .59194 m .35636 .65466 L .5 .74578 L closepath p
.029 .388 .849 r
F P
s
P
p
.002 w
.4168 .39477 m .45578 .29834 L .3433 .351 L closepath p
.168 0 .005 r
F P
s
P
p
.002 w
.30262 .45997 m .4168 .39477 L .3433 .351 L closepath p
.057 0 .209 r
F P
s
P
p
.002 w
.74854 .4729 m .61499 .42179 L .59728 .5671 L closepath p
.912 .736 .65 r
F P
s
P
p
.002 w
.74854 .4729 m .61499 .42179 L .59005 .31539 L closepath p
.89 .508 .307 r
F P
s
P
p
.002 w
.65265 .63107 m .74854 .4729 L .59728 .5671 L closepath p
.937 .888 .758 r
F P
s
P
p
.002 w
.45578 .29834 m .4168 .39477 L .59005 .31539 L closepath p
.588 .142 .159 r
F P
s
P
p
.002 w
.59728 .5671 m .65265 .63107 L .5 .74578 L closepath p
.822 .894 .893 r
F P
s
P
p
.002 w
.45152 .59845 m .5 .74578 L .33138 .59194 L closepath p
.452 .58 .891 r
F P
s
P
p
.002 w
.45152 .59845 m .59728 .5671 L .5 .74578 L closepath p
.717 .749 .877 r
F P
s
P
p
.002 w
.33138 .59194 m .30262 .45997 L .4168 .39477 L closepath p
.28 .235 .628 r
F P
s
P
p
.002 w
.59005 .31539 m .61499 .42179 L .4168 .39477 L closepath p
.772 .428 .404 r
F P
s
P
p
.002 w
.45152 .59845 m .4168 .39477 L .33138 .59194 L closepath p
.493 .477 .768 r
F P
s
P
p
.002 w
.59728 .5671 m .61499 .42179 L .4168 .39477 L closepath p
.797 .625 .668 r
F P
s
P
p
.002 w
.45152 .59845 m .59728 .5671 L .4168 .39477 L closepath p
.715 .633 .765 r
F P
s
P
p
P
p
P
% End of Graphics
MathPictureEnd

:[font = subsection; inactive; Cclosed; preserveAspect; startGroup]
Draw a  Polytope  by a New VisibleSkeleton Graphics
:[font = input; preserveAspect; startGroup; animationSpeed = 80; infiniteLoop; loopDistance = 1]
Do[ uvp=circle[a,5,2];
	Show[skel3D[uvp], 
		Boxed->False, 
 		ViewPoint -> getMmaViewPoint[uvp,skel3D[uvp]],
 		SphericalRegion->True
 	],
 	{a,Pi/11,Pi,Pi/18}
 ]
:[font = postscript; PostScript; formatAsPostScript; output; inactive; preserveAspect; pictureLeft = 34; pictureWidth = 282; pictureHeight = 282; animationSpeed = 6]
%!
%%Creator: Mathematica
%%AspectRatio: 1 
MathPictureStart
%% Graphics3D
/Courier findfont 10  scalefont  setfont
% Scaling calculations
0.02381 0.952381 0.02381 0.952381 [
[ 0 0 0 0 ]
[ 1 1 0 0 ]
] MathScale
% Start of Graphics
1 setlinecap
1 setlinejoin
newpath
[ ] 0 setdash
0 g
0 0 m
1 0 L
1 1 L
0 1 L
closepath
clip
newpath
p
P
p
[ .01 .012 ] 0 setdash
.0035 w
.55958 .57536 m
.42213 .65964 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.41118 .56041 m
.55958 .57536 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.63995 .42369 m
.55958 .57536 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.53652 .42583 m
.55958 .57536 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.67152 .53478 m
.55958 .57536 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.42683 .44953 m
.55958 .57536 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.64938 .64203 m
.55958 .57536 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.53938 .67959 m
.55958 .57536 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.42683 .44953 m
.36746 .3862 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.42683 .44953 m
.27448 .52459 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.42683 .44953 m
.41118 .56041 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.42683 .44953 m
.53652 .42583 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.5 .28594 m
.42683 .44953 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.53652 .42583 m
.63995 .42369 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.5 .28594 m
.53652 .42583 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.41118 .56041 m
.27448 .52459 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.41118 .56041 m
.42213 .65964 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.63995 .42369 m
.64017 .34908 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.74854 .4729 m
.63995 .42369 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.67152 .53478 m
.63995 .42369 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.5 .28594 m
.63995 .42369 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.42213 .65964 m
.27448 .52459 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.35636 .65466 m
.42213 .65964 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.5 .74578 m
.42213 .65964 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.53938 .67959 m
.42213 .65964 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.27448 .52459 m
.36746 .3862 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.3433 .351 m
.36746 .3862 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.5 .28594 m
.36746 .3862 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.67152 .53478 m
.74854 .4729 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.64938 .64203 m
.67152 .53478 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.5 .28594 m
.3433 .351 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.5 .28594 m
.64017 .34908 L
s
P
p
.0085 w
.5 .28594 m
.45578 .29834 L
s
P
p
.0085 w
.5 .28594 m
.59005 .31539 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.53938 .67959 m
.64938 .64203 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.53938 .67959 m
.5 .74578 L
s
P
p
.0085 w
.3433 .351 m
.27448 .52459 L
s
P
p
.0085 w
.33138 .59194 m
.27448 .52459 L
s
P
p
.0085 w
.35636 .65466 m
.27448 .52459 L
s
P
p
.0085 w
.30262 .45997 m
.27448 .52459 L
s
P
p
.0085 w
.64938 .64203 m
.74854 .4729 L
s
P
p
.0085 w
.64938 .64203 m
.65265 .63107 L
s
P
p
.0085 w
.5 .74578 m
.64938 .64203 L
s
P
p
.0085 w
.74854 .4729 m
.64017 .34908 L
s
P
p
.0085 w
.59005 .31539 m
.64017 .34908 L
s
P
p
.0085 w
.35636 .65466 m
.33138 .59194 L
s
P
p
.0085 w
.5 .74578 m
.35636 .65466 L
s
P
p
.0085 w
.4168 .39477 m
.3433 .351 L
s
P
p
.0085 w
.45578 .29834 m
.3433 .351 L
s
P
p
.0085 w
.30262 .45997 m
.3433 .351 L
s
P
p
.0085 w
.74854 .4729 m
.59728 .5671 L
s
P
p
.0085 w
.61499 .42179 m
.74854 .4729 L
s
P
p
.0085 w
.65265 .63107 m
.74854 .4729 L
s
P
p
.0085 w
.59005 .31539 m
.74854 .4729 L
s
P
p
.0085 w
.45578 .29834 m
.4168 .39477 L
s
P
p
.0085 w
.59005 .31539 m
.45578 .29834 L
s
P
p
.0085 w
.5 .74578 m
.45152 .59845 L
s
P
p
.0085 w
.5 .74578 m
.59728 .5671 L
s
P
p
.0085 w
.5 .74578 m
.33138 .59194 L
s
P
p
.0085 w
.5 .74578 m
.65265 .63107 L
s
P
p
.0085 w
.30262 .45997 m
.4168 .39477 L
s
P
p
.0085 w
.30262 .45997 m
.33138 .59194 L
s
P
p
.0085 w
.65265 .63107 m
.59728 .5671 L
s
P
p
.0085 w
.59005 .31539 m
.4168 .39477 L
s
P
p
.0085 w
.59005 .31539 m
.61499 .42179 L
s
P
p
.0085 w
.33138 .59194 m
.4168 .39477 L
s
P
p
.0085 w
.33138 .59194 m
.45152 .59845 L
s
P
p
.0085 w
.61499 .42179 m
.4168 .39477 L
s
P
p
.0085 w
.61499 .42179 m
.59728 .5671 L
s
P
p
.0085 w
.45152 .59845 m
.4168 .39477 L
s
P
p
.0085 w
.59728 .5671 m
.45152 .59845 L
s
P
p
.0085 w
.59728 .5671 m
.4168 .39477 L
s
P
p
P
p
P
% End of Graphics
MathPictureEnd

:[font = postscript; PostScript; formatAsPostScript; output; inactive; preserveAspect; pictureLeft = 34; pictureWidth = 282; pictureHeight = 282; endGroup; endGroup; endGroup; animationSpeed = 6]
%!
%%Creator: Mathematica
%%AspectRatio: 1 
MathPictureStart
%% Graphics3D
/Courier findfont 10  scalefont  setfont
% Scaling calculations
0.02381 0.952381 0.02381 0.952381 [
[ 0 0 0 0 ]
[ 1 1 0 0 ]
] MathScale
% Start of Graphics
1 setlinecap
1 setlinejoin
newpath
[ ] 0 setdash
0 g
0 0 m
1 0 L
1 1 L
0 1 L
closepath
clip
newpath
p
P
p
[ .01 .012 ] 0 setdash
.0035 w
.5948 .57104 m
.44849 .66257 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.44131 .56463 m
.5948 .57104 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.65383 .41268 m
.5948 .57104 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.5579 .42232 m
.5948 .57104 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.68931 .52371 m
.5948 .57104 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.45157 .45384 m
.5948 .57104 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.6644 .63449 m
.5948 .57104 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.56246 .6774 m
.5948 .57104 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.45157 .45384 m
.37815 .39616 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.45157 .45384 m
.29488 .53776 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.45157 .45384 m
.44131 .56463 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.45157 .45384 m
.5579 .42232 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.5 .28594 m
.45157 .45384 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.44131 .56463 m
.29488 .53776 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.44131 .56463 m
.44849 .66257 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.5579 .42232 m
.65383 .41268 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.5 .28594 m
.5579 .42232 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.29488 .53776 m
.37815 .39616 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.33413 .3643 m
.37815 .39616 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.5 .28594 m
.37815 .39616 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.44849 .66257 m
.29488 .53776 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.36811 .66095 m
.44849 .66257 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.5 .74578 m
.44849 .66257 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.56246 .6774 m
.44849 .66257 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.33413 .3643 m
.29488 .53776 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.36811 .66095 m
.29488 .53776 L
s
P
p
.0085 w
.32178 .60116 m
.29488 .53776 L
s
P
p
.0085 w
.29212 .47396 m
.29488 .53776 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.65383 .41268 m
.6331 .33758 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.73913 .45598 m
.65383 .41268 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.68931 .52371 m
.65383 .41268 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.5 .28594 m
.65383 .41268 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.5 .28594 m
.33413 .3643 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.5 .28594 m
.6331 .33758 L
s
P
p
.0085 w
.5 .28594 m
.43078 .30341 L
s
P
p
.0085 w
.5 .28594 m
.5606 .30876 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.56246 .6774 m
.6644 .63449 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.56246 .6774 m
.5 .74578 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.68931 .52371 m
.73913 .45598 L
s
P
p
[ .01 .012 ] 0 setdash
.0035 w
.6644 .63449 m
.68931 .52371 L
s
P
p
.0085 w
.36811 .66095 m
.32178 .60116 L
s
P
p
.0085 w
.5 .74578 m
.36811 .66095 L
s
P
p
.0085 w
.6644 .63449 m
.73913 .45598 L
s
P
p
.0085 w
.6644 .63449 m
.64517 .62371 L
s
P
p
.0085 w
.5 .74578 m
.6644 .63449 L
s
P
p
.0085 w
.37056 .40301 m
.33413 .3643 L
s
P
p
.0085 w
.43078 .30341 m
.33413 .3643 L
s
P
p
.0085 w
.29212 .47396 m
.33413 .3643 L
s
P
p
.0085 w
.73913 .45598 m
.6331 .33758 L
s
P
p
.0085 w
.5606 .30876 m
.6331 .33758 L
s
P
p
.0085 w
.29212 .47396 m
.37056 .40301 L
s
P
p
.0085 w
.29212 .47396 m
.32178 .60116 L
s
P
p
.0085 w
.43078 .30341 m
.37056 .40301 L
s
P
p
.0085 w
.5606 .30876 m
.43078 .30341 L
s
P
p
.0085 w
.5 .74578 m
.42417 .60173 L
s
P
p
.0085 w
.5 .74578 m
.56552 .56245 L
s
P
p
.0085 w
.5 .74578 m
.32178 .60116 L
s
P
p
.0085 w
.5 .74578 m
.64517 .62371 L
s
P
p
.0085 w
.32178 .60116 m
.37056 .40301 L
s
P
p
.0085 w
.32178 .60116 m
.42417 .60173 L
s
P
p
.0085 w
.73913 .45598 m
.56552 .56245 L
s
P
p
.0085 w
.57758 .41456 m
.73913 .45598 L
s
P
p
.0085 w
.64517 .62371 m
.73913 .45598 L
s
P
p
.0085 w
.5606 .30876 m
.73913 .45598 L
s
P
p
.0085 w
.5606 .30876 m
.37056 .40301 L
s
P
p
.0085 w
.5606 .30876 m
.57758 .41456 L
s
P
p
.0085 w
.64517 .62371 m
.56552 .56245 L
s
P
p
.0085 w
.42417 .60173 m
.37056 .40301 L
s
P
p
.0085 w
.56552 .56245 m
.42417 .60173 L
s
P
p
.0085 w
.57758 .41456 m
.37056 .40301 L
s
P
p
.0085 w
.57758 .41456 m
.56552 .56245 L
s
P
p
.0085 w
.56552 .56245 m
.37056 .40301 L
s
P
p
P
p
P
% End of Graphics
MathPictureEnd

:[font = section; inactive; Cclosed; preserveAspect; startGroup]
Unfolding the Polytope
:[font = input; preserveAspect; endGroup; animationSpeed = 23; infiniteLoop; loopDistance = 1]
UnfoldPolytope[facets1] 
:[font = section; inactive; Cclosed; preserveAspect; startGroup]
Disconnect cddmathlink
:[font = input; preserveAspect; startGroup]
Uninstall[cddml]
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
Uninstall[$Failed]
;[o]
Uninstall[$Failed]
^*)
