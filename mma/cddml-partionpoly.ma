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
	paletteColors = 128; currentKernel; 
]
:[font = title; inactive; preserveAspect]
cddmathlink
Onn's Partition Polytope by MathLink to cddlib
by Komei Fukuda
April 13, 2000
;[s]
6:0,0;11,1;40,2;48,3;52,4;58,5;89,-1;
6:1,21,16,Times,1,24,3651,9032,33889;1,21,16,Times,1,24,0,0,0;1,22,17,Times,3,24,960,30237,6342;1,21,16,Times,1,24,0,0,0;1,21,16,Times,1,24,33889,1793,1793;1,21,16,Times,1,24,0,0,0;
:[font = subsection; inactive; initialization; preserveAspect; startGroup]
Connecting  cddmathlink
:[font = text; inactive; initialization; preserveAspect; startGroup]
You just put the compiled cddmathlink for your computer in some directory.  In this example, the name of the directory is "~/Math".
;[s]
2:0,0;122,1;131,-1;
2:1,11,8,Times,0,12,0,0,0;1,10,8,Courier,1,12,0,0,0;
:[font = input; initialization; preserveAspect]
*)
Off[General::spell1]; Off[General::spell];
(*
:[font = input; initialization; preserveAspect; startGroup]
*)
cddml=
Install["~/Math/cddmathlink"]
(*
:[font = output; output; inactive; initialization; preserveAspect; endGroup; endGroup; endGroup]
LinkObject["'/Users/fukuda/Math/cddmathlink'", 2, 2]
;[o]
LinkObject['/Users/fukuda/Math/cddmathlink', 2, 2]
:[font = subsection; inactive; preserveAspect]
Partition Polytope 
:[font = text; inactive; preserveAspect]
For any fixed positive integer p, d and n (p<n), fix a d times n matrix A.
:[font = input; preserveAspect]
<<~/Math/Combinatorica3.m
:[font = input; preserveAspect; startGroup]
p=3; d=2; n=3;  
MatrixForm[A=Table[Table[Random[Integer,{0,3}],{n}],{d}]]
:[font = output; output; inactive; preserveAspect; endGroup]
MatrixForm[{{2, 1, 0}, {2, 1, 2}}]
;[o]
2   1   0

2   1   2
:[font = input; preserveAspect; startGroup]
Show[Graphics[Map[Point[Transpose[A][[#]]]&, Range[n]]]]
:[font = postscript; PostScript; formatAsPostScript; output; inactive; preserveAspect; pictureLeft = 34; pictureWidth = 282; pictureHeight = 174]
%!
%%Creator: Mathematica
%%AspectRatio: .61803 
MathPictureStart
/Mabs {
Mgmatrix idtransform
Mtmatrix dtransform
} bind def
/Mabsadd { Mabs
3 -1 roll add
3 1 roll add
exch } bind def
%% Graphics
%%IncludeResource: font Courier
%%IncludeFont: Courier
/Courier findfont 10  scalefont  setfont
% Scaling calculations
0.02381 0.47619 -0.573889 0.588604 [
[ 0 0 0 0 ]
[ 1 .61803 0 0 ]
] MathScale
% Start of Graphics
1 setlinecap
1 setlinejoin
newpath
0 0 m
1 0 L
1 .61803 L
0 .61803 L
closepath
clip
newpath
0 g
.008 w
.97619 .60332 Mdot
.5 .01472 Mdot
.02381 .60332 Mdot
% End of Graphics
MathPictureEnd

:[font = output; output; inactive; preserveAspect; endGroup]
Graphics["<<>>"]
;[o]
-Graphics-
:[font = input; preserveAspect; startGroup]
s={};
Do[
	Do[AppendTo[s,{i,j}],
		{i,1,p}],
	{j,1,n}
];
s

:[font = output; output; inactive; preserveAspect; endGroup]
{{1, 1}, {2, 1}, {3, 1}, {1, 2}, {2, 2}, {3, 2}, {1, 3}, 
 
  {2, 3}, {3, 3}}
;[o]
{{1, 1}, {2, 1}, {3, 1}, {1, 2}, {2, 2}, {3, 2}, {1, 3}, 
 
  {2, 3}, {3, 3}}
:[font = subsection; inactive; preserveAspect; startGroup]
Generating the Graph Structure
:[font = input; preserveAspect; startGroup]
?AllVerticesWithAdjacency
:[font = print; inactive; preserveAspect]
AllVerticesWithAdjacency[m,d+1,A] generates all extreme points
   (vertices) and extreme rays of the convex polyhedron in R^(d+1) given
   as the solution set to an inequality system  A x >= 0 where   A is an
   m*(d+1) matrix  and x=(1,x1,...,xd). The output is {{extlist,
   linearity}, inclist, adjlist} where extlist , inclist, adjlist are
   the extreme point list, the incidence list, the adjacency list (of
   extreme points and rays). Each vertex (ray) has the first component 1
   (0). If the convex polyhedron has a vertex, linearity is the empty
   set.  If the convex polyhedron has no vertices, extlist is a
   (nonunique) set of generators of the polyhedron where those
   generators in the linearity list are considered as linearity space
   (of points satisfying A (0, x1, x2, ...., xd) = 0) generators.
:[font = input; preserveAspect; startGroup]
{{vertices,linearity},incidences,adjacencies}=
	AllVerticesWithAdjacency[m,d1,Flatten[a]]
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
{{{{1., 0, 0, 1.}, {1., 0, 1., 1.}, {1., 1., 1., 1.}, {1., 1., 0, 1.}, 
 
    {1., 1., 1., 0}, {1., 1., 0, 0}, {1., 0, 1., 0}, {1., 0, 0, 0}}, 
 
   {}}, {{1, 2, 6}, {1, 5, 6}, {4, 5, 6}, {2, 4, 6}, {3, 4, 5}, 
 
   {2, 3, 4}, {1, 3, 5}, {1, 2, 3}}, 
 
  {{2, 4, 8}, {1, 3, 7}, {2, 4, 5}, {1, 3, 6}, {3, 6, 7}, {4, 5, 8}, 
 
   {2, 5, 8}, {1, 6, 7}}}
;[o]
{{{{1., 0, 0, 1.}, {1., 0, 1., 1.}, {1., 1., 1., 1.}, {1., 1., 0, 1.}, 
 
    {1., 1., 1., 0}, {1., 1., 0, 0}, {1., 0, 1., 0}, {1., 0, 0, 0}}, 
 
   {}}, {{1, 2, 6}, {1, 5, 6}, {4, 5, 6}, {2, 4, 6}, {3, 4, 5}, 
 
   {2, 3, 4}, {1, 3, 5}, {1, 2, 3}}, 
 
  {{2, 4, 8}, {1, 3, 7}, {2, 4, 5}, {1, 3, 6}, {3, 6, 7}, {4, 5, 8}, 
 
   {2, 5, 8}, {1, 6, 7}}}
:[font = text; inactive; preserveAspect; endGroup]
The graph structure is output as the adjacency  list  adjacencies. For example, the first list {2, 4 ,8} represents the neighbour vertices of the first vertex 1.
;[s]
3:0,0;54,1;67,2;163,-1;
3:1,11,8,Times,0,12,0,0,0;1,10,8,Courier,1,12,0,0,0;1,11,8,Times,0,12,0,0,0;
:[font = subsection; inactive; Cclosed; preserveAspect; startGroup]
Convex Hull (Facet Generation)
:[font = input; preserveAspect; startGroup]
?AllFacets
:[font = print; inactive; preserveAspect]
AllFacets[n,d+1,G] generates all facet inequalities of the convex
   polyhedron in R^(d+1) generated by points and rays given in the rows
   of an n*(d+1) matrix G.  Each point (ray) must have 1 (0) in the
   first coordinate.  The output is {{faclist, equalities}, inclist}
   where faclist is  the facet  list and inclist is the incidence list. 
   If the convex polyhedron is not full-dimensional, extlist is a
   (nonunique) set of inequalities of the polyhedron where those
   inequalities in the equalities list are considered as equalities.
:[font = text; inactive; preserveAspect]
We have computed all the vertices of a 3-cube.  Let's try the reverse operation.  First check the size of the list of vertics.   It should reconstruct the facets we have started with.
:[font = input; preserveAspect; startGroup]
{n, d1}=Dimensions[vertices]
:[font = output; output; inactive; preserveAspect; endGroup]
{8, 4}
;[o]
{8, 4}
:[font = input; preserveAspect; startGroup]

{{facets,equalities}, fincidences}= AllFacets[n,d1,Flatten[vertices]]
:[font = output; output; inactive; preserveAspect; endGroup]
{{{{0, 0, 1., 0}, {0, 0, 0, 1.}, {0, 1., 0, 0}, {1., 0, -1., 0}, 
 
    {1., -1., 0, 0}, {1., 0, 0, -1.}}, {}}, 
 
  {{1, 4, 6, 8}, {5, 6, 7, 8}, {1, 2, 7, 8}, {2, 3, 5, 7}, 
 
   {3, 4, 5, 6}, {1, 2, 3, 4}}}
;[o]
{{{{0, 0, 1., 0}, {0, 0, 0, 1.}, {0, 1., 0, 0}, {1., 0, -1., 0}, 
 
    {1., -1., 0, 0}, {1., 0, 0, -1.}}, {}}, 
 
  {{1, 4, 6, 8}, {5, 6, 7, 8}, {1, 2, 7, 8}, {2, 3, 5, 7}, 
 
   {3, 4, 5, 6}, {1, 2, 3, 4}}}
:[font = text; inactive; preserveAspect]
We can compute how the facets are connected by using AllFacetsWithAdjacency function.
;[s]
3:0,0;53,1;75,2;87,-1;
3:1,11,8,Times,0,12,0,0,0;1,10,8,Times,1,12,0,0,0;1,11,8,Times,0,12,0,0,0;
:[font = input; preserveAspect; startGroup]

{{facets,equalities}, fincidences,fadjacencies}= 
	AllFacetsWithAdjacency[n,d1,Flatten[vertices]]
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
{{{{0, 0, 1., 0}, {0, 0, 0, 1.}, {0, 1., 0, 0}, {1., 0, -1., 0}, 
 
    {1., -1., 0, 0}, {1., 0, 0, -1.}}, {}}, 
 
  {{1, 4, 6, 8}, {5, 6, 7, 8}, {1, 2, 7, 8}, {2, 3, 5, 7}, 
 
   {3, 4, 5, 6}, {1, 2, 3, 4}}, 
 
  {{2, 3, 5, 6}, {1, 3, 4, 5}, {1, 2, 4, 6}, {2, 3, 5, 6}, 
 
   {1, 2, 4, 6}, {1, 3, 4, 5}}}
;[o]
{{{{0, 0, 1., 0}, {0, 0, 0, 1.}, {0, 1., 0, 0}, {1., 0, -1., 0}, 
 
    {1., -1., 0, 0}, {1., 0, 0, -1.}}, {}}, 
 
  {{1, 4, 6, 8}, {5, 6, 7, 8}, {1, 2, 7, 8}, {2, 3, 5, 7}, 
 
   {3, 4, 5, 6}, {1, 2, 3, 4}}, 
 
  {{2, 3, 5, 6}, {1, 3, 4, 5}, {1, 2, 4, 6}, {2, 3, 5, 6}, 
 
   {1, 2, 4, 6}, {1, 3, 4, 5}}}
:[font = text; inactive; preserveAspect; startGroup]
If you want to compute an inequality description of the one-dimensional cone in R^3 with a vertex at origin and containing the direction (1,1,1), you must set up the input (generator) data as:
:[font = input; preserveAspect; startGroup]
coneGenerators={{1,0,0,0},{0,1,1,1}}
:[font = output; output; inactive; preserveAspect; endGroup]
{{1, 0, 0, 0}, {0, 1, 1, 1}}
;[o]
{{1, 0, 0, 0}, {0, 1, 1, 1}}
:[font = input; preserveAspect; startGroup]
{{cfacets,cequalities}, cfincidences}= 
	AllFacets[2,4,Flatten[coneGenerators]]
:[font = output; output; inactive; preserveAspect; endGroup]
{{{{1., 0, 0, 0}, {0, 1., 0, 0}, {0, -1., 1., 0}, {0, -1., 0, 1.}}, 
 
   {3, 4}}, {{2}, {1}}}
;[o]
{{{{1., 0, 0, 0}, {0, 1., 0, 0}, {0, -1., 1., 0}, {0, -1., 0, 1.}}, 
 
   {3, 4}}, {{2}, {1}}}
:[font = text; inactive; preserveAspect; endGroup; endGroup]
Since the equalities list contains 3 and 4, of the four output inequalities , the third and the forth must be considered as equalities.  It is important to note that this cone can have infinitely many different minimal inequality descriptions, since it is not full-dimensional.
:[font = subsection; inactive; Cclosed; preserveAspect; startGroup]
A Larger Example (Random 0/1 Polytopes)
:[font = text; inactive; preserveAspect; fontSize = 13; fontName = "Times"]
Let's compute the convex hull of 0/1 points in R^d.  First generate 0/1 points.  Below each point with the first component 0 is considered as a direction that must be included in the convex hull.
:[font = input; preserveAspect]
n=30; d1=8;
:[font = input; preserveAspect; startGroup]
points=Table[Prepend[Table[Random[Integer,{0,1}],{d1-1}],0],{n}]
:[font = output; output; inactive; preserveAspect; endGroup]
{{0, 1, 0, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 1, 1}, 
 
  {0, 0, 1, 1, 0, 0, 1, 0}, {0, 1, 0, 0, 1, 1, 0, 0}, 
 
  {0, 1, 1, 1, 0, 1, 1, 1}, {0, 0, 1, 0, 1, 1, 1, 1}, 
 
  {0, 0, 0, 1, 0, 0, 0, 0}, {0, 0, 0, 1, 1, 0, 1, 0}, 
 
  {0, 1, 0, 0, 0, 1, 0, 0}, {0, 1, 1, 0, 1, 1, 0, 1}, 
 
  {0, 1, 1, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0}, 
 
  {0, 0, 1, 1, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 1, 1, 0}, 
 
  {0, 1, 0, 1, 1, 1, 1, 1}, {0, 0, 1, 0, 0, 0, 0, 0}, 
 
  {0, 0, 0, 1, 1, 0, 1, 0}, {0, 0, 0, 1, 1, 0, 1, 1}, 
 
  {0, 1, 1, 0, 0, 0, 0, 1}, {0, 0, 1, 0, 1, 1, 1, 0}, 
 
  {0, 0, 0, 1, 0, 1, 1, 1}, {0, 1, 0, 0, 0, 0, 1, 0}, 
 
  {0, 0, 1, 0, 0, 1, 1, 0}, {0, 0, 0, 1, 1, 0, 0, 1}, 
 
  {0, 1, 1, 1, 0, 1, 1, 0}, {0, 0, 0, 1, 0, 0, 0, 1}, 
 
  {0, 1, 1, 0, 0, 1, 1, 0}, {0, 1, 1, 0, 1, 1, 0, 1}, 
 
  {0, 0, 1, 1, 0, 0, 0, 1}, {0, 1, 1, 1, 1, 1, 0, 0}}
;[o]
{{0, 1, 0, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 1, 1}, 
 
  {0, 0, 1, 1, 0, 0, 1, 0}, {0, 1, 0, 0, 1, 1, 0, 0}, 
 
  {0, 1, 1, 1, 0, 1, 1, 1}, {0, 0, 1, 0, 1, 1, 1, 1}, 
 
  {0, 0, 0, 1, 0, 0, 0, 0}, {0, 0, 0, 1, 1, 0, 1, 0}, 
 
  {0, 1, 0, 0, 0, 1, 0, 0}, {0, 1, 1, 0, 1, 1, 0, 1}, 
 
  {0, 1, 1, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0}, 
 
  {0, 0, 1, 1, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 1, 1, 0}, 
 
  {0, 1, 0, 1, 1, 1, 1, 1}, {0, 0, 1, 0, 0, 0, 0, 0}, 
 
  {0, 0, 0, 1, 1, 0, 1, 0}, {0, 0, 0, 1, 1, 0, 1, 1}, 
 
  {0, 1, 1, 0, 0, 0, 0, 1}, {0, 0, 1, 0, 1, 1, 1, 0}, 
 
  {0, 0, 0, 1, 0, 1, 1, 1}, {0, 1, 0, 0, 0, 0, 1, 0}, 
 
  {0, 0, 1, 0, 0, 1, 1, 0}, {0, 0, 0, 1, 1, 0, 0, 1}, 
 
  {0, 1, 1, 1, 0, 1, 1, 0}, {0, 0, 0, 1, 0, 0, 0, 1}, 
 
  {0, 1, 1, 0, 0, 1, 1, 0}, {0, 1, 1, 0, 1, 1, 0, 1}, 
 
  {0, 0, 1, 1, 0, 0, 0, 1}, {0, 1, 1, 1, 1, 1, 0, 0}}
:[font = input; preserveAspect; startGroup]
Dimensions[points]
:[font = output; output; inactive; preserveAspect; endGroup]
{30, 8}
;[o]
{30, 8}
:[font = input; preserveAspect]

{CPUtime, {{facets,equalities}, inc, adj}}= 
	Timing[AllFacetsWithAdjacency[n,d1,Flatten[points]]];
:[font = input; preserveAspect; startGroup]
{CPUtime,Length[facets]}
:[font = output; output; inactive; preserveAspect; endGroup]
{0.01666666666666661*Second, 31}
;[o]
{0.0166667 Second, 31}
:[font = text; inactive; preserveAspect]
Usually facets of 0/1 polytopes are very pretty and their coefficients are small integers.
:[font = input; preserveAspect; startGroup]
Take[facets,5]
:[font = output; output; inactive; preserveAspect; endGroup]
{{0, 0, 0, 1., -1., 1., 0, 0}, {0, 2., 0, 1., -1., -1., 2., 0}, 
 
  {0, 2., 0, 1., 0, -1., 2., -1.}, {0, 1., 0, 1., 0, 0, 1., -1.}, 
 
  {0, 1., 1., 0, 1., 1., -1., 0}}
;[o]
{{0, 0, 0, 1., -1., 1., 0, 0}, {0, 2., 0, 1., -1., -1., 2., 0}, 
 
  {0, 2., 0, 1., 0, -1., 2., -1.}, {0, 1., 0, 1., 0, 0, 1., -1.}, 
 
  {0, 1., 1., 0, 1., 1., -1., 0}}
:[font = text; inactive; preserveAspect]
Also each facet tends to have many neighbours.
:[font = input; preserveAspect; startGroup]
Take[adj,5]
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
{{2, 4, 6, 7, 9, 10, 11, 14, 21, 23, 26, 30, 31}, 
 
  {1, 3, 10, 11, 26, 28, 29}, {2, 4, 10, 26, 27, 28}, 
 
  {1, 3, 6, 10, 26, 27, 31}, {6, 9, 10, 12, 13, 21, 24, 31}}
;[o]
{{2, 4, 6, 7, 9, 10, 11, 14, 21, 23, 26, 30, 31}, 
 
  {1, 3, 10, 11, 26, 28, 29}, {2, 4, 10, 26, 27, 28}, 
 
  {1, 3, 6, 10, 26, 27, 31}, {6, 9, 10, 12, 13, 21, 24, 31}}
:[font = subsection; inactive; Cclosed; preserveAspect; startGroup]
Disconnecting  cddmathlink
:[font = input; preserveAspect; startGroup]
Uninstall[mlcdd]
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
"'/Users/fukuda/Math/cddmathlink'"
;[o]
'/Users/fukuda/Math/cddmathlink'
^*)
