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
:[font = title; inactive; preserveAspect; fontSize = 27; fontName = "Times"]
Diet Problem
An Application of Vertex Enumeration
with MathLink to cddlib
;[s]
5:0,0;12,1;55,2;63,3;67,4;73,-1;
5:1,24,18,Times,1,27,3651,9032,33889;1,24,18,Times,1,27,0,0,0;1,22,17,Times,3,24,960,30237,6342;1,24,18,Times,1,27,0,0,0;1,24,18,Times,1,27,33889,1793,1793;
:[font = subtitle; inactive; preserveAspect]
Komei Fukuda, fukuda@ifor.math.ethz.ch
Swiss Federal Institute of Technology, Lausanne and Zurich
March 14, 1999
;[s]
1:0,0;113,-1;
1:1,15,12,Times,3,17,0,0,0;
:[font = section; inactive; initialization; Cclosed; preserveAspect; startGroup]
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
cddml=Install["~/Math/cddmathlink"]
(*
:[font = output; output; inactive; initialization; preserveAspect; endGroup; endGroup; endGroup]
LinkObject["~/Math/cddmathlink", 21, 21]
;[o]
LinkObject[~/Math/cddmathlink, 21, 21]
:[font = section; inactive; Cclosed; preserveAspect; fontSize = 20; fontName = "Times"; startGroup]
What is Diet Problem?
:[font = text; inactive; preserveAspect]
The following diet problem is taken from V. Chvatal's  great book on Linear Programming ("Linear Programming", W.H.Freeman and Company,1983).   It is to design a cheapest meal with six possible items below to satisfy prescribed nutritional needs.  Please see Page 3 of the book.
:[font = input; preserveAspect; fontSize = 14; fontName = "Courier"; startGroup]
var={"","Oatmeal","Chicken","Eggs","Milk","Cherry Pie", 
	"Pork Beans"};

price={"Price/Ser", "3c", "24c", "13c", "9c", "20c", "19c"}
:[font = output; output; inactive; preserveAspect; endGroup]
{"Price/Ser", "3c", "24c", "13c", "9c", "20c", "19c"}
;[o]
{Price/Ser, 3c, 24c, 13c, 9c, 20c, 19c}
:[font = input; preserveAspect; fontSize = 14; fontName = "Courier"]
MatrixForm[dietproblem1=
{{0, 1, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0},
 {0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0}, 
 {0, 0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 0, 1}, 
 {4, -1, 0, 0, 0, 0, 0}, {3, 0, -1, 0, 0, 0, 0},
 {2, 0, 0, -1, 0, 0, 0}, {8, 0, 0, 0, -1, 0, 0}, 
 {2, 0, 0, 0, 0, -1, 0}, {2, 0, 0, 0, 0, 0, -1},
 {-2000, 110, 205, 160, 160, 420, 260}, 
 {-55, 4, 32, 13, 8, 4, 14}, 
 {-800, 2, 12, 54, 285, 22, 80}}];
:[font = input; preserveAspect; startGroup]
TableForm[table1=Prepend[Prepend[dietproblem1,var],price]]

:[font = output; output; inactive; preserveAspect; endGroup]
TableForm[{{"Price/Ser", "3c", "24c", "13c", "9c", "20c", "19c"}, 
 
   {"", "Oatmeal", "Chicken", "Eggs", "Milk", "Cherry Pie", "Pork Beans"}, 
 
   {0, 1, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0}, 
 
   {0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 0, 1}, 
 
   {4, -1, 0, 0, 0, 0, 0}, {3, 0, -1, 0, 0, 0, 0}, {2, 0, 0, -1, 0, 0, 0}, 
 
   {8, 0, 0, 0, -1, 0, 0}, {2, 0, 0, 0, 0, -1, 0}, {2, 0, 0, 0, 0, 0, -1}, 
 
   {-2000, 110, 205, 160, 160, 420, 260}, {-55, 4, 32, 13, 8, 4, 14}, 
 
   {-800, 2, 12, 54, 285, 22, 80}}]
;[o]
Price/Ser   3c        24c       13c    9c     20c          19c

            Oatmeal   Chicken   Eggs   Milk   Cherry Pie   Pork Beans

0           1         0         0      0      0            0

0           0         1         0      0      0            0

0           0         0         1      0      0            0

0           0         0         0      1      0            0

0           0         0         0      0      1            0

0           0         0         0      0      0            1

4           -1        0         0      0      0            0

3           0         -1        0      0      0            0

2           0         0         -1     0      0            0

8           0         0         0      -1     0            0

2           0         0         0      0      -1           0

2           0         0         0      0      0            -1

-2000       110       205       160    160    420          260

-55         4         32        13     8      4            14

-800        2         12        54     285    22           80
:[font = input; preserveAspect; fontSize = 14; fontName = "Courier"]
m=Transpose[Drop[Transpose[dietproblem1],1]];
b=-First[Transpose[dietproblem1]];
c={3, 24, 13, 9, 20, 19};
:[font = text; inactive; preserveAspect]
By using the build-in LP optimizer of Mathematica, one can easily compute the optimal solution.
;[s]
3:0,0;38,1;49,2;95,-1;
3:1,11,8,Times,0,12,0,0,0;1,10,8,Times,2,12,0,0,0;1,11,8,Times,0,12,0,0,0;
:[font = input; preserveAspect; fontSize = 14; fontName = "Courier"; startGroup]
lps=LinearProgramming[c, m,b]
:[font = output; output; inactive; preserveAspect; fontSize = 14; fontName = "Courier"; endGroup]
{4, 0, 0, 9/2, 2, 0}
;[o]
          9
{4, 0, 0, -, 2, 0}
          2
:[font = input; preserveAspect; fontSize = 14; fontName = "Courier"; startGroup]
optvalue= N[c.lps]
:[font = output; output; inactive; preserveAspect; fontSize = 14; fontName = "Courier"; endGroup]
92.5
;[o]
92.5
:[font = text; inactive; preserveAspect]
We can see the optimal solution better in the following table.   It is certainly not an exciting menu.   In fact, an optimal solution to any optimization problem tends to be extreme, and thus it must be modified for practical purposes.
:[font = input; preserveAspect; fontSize = 14; fontName = "Courier"; startGroup]
TableForm[Join[{var},{Prepend[N[lps],optvalue]}]]

:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
TableForm[{{"", "Oatmeal", "Chicken", "Eggs", "Milk", "Cherry Pie", "Pork\
 
     Beans"}, {92.5, 4., 0, 0, 4.5, 2., 0}}]
;[o]
       Oatmeal   Chicken   Eggs   Milk   Cherry Pie   Pork Beans

92.5   4.        0         0      4.5    2.           0
:[font = section; inactive; Cclosed; preserveAspect; startGroup]
Why is the Vertex Enumeration Useful?
:[font = subsection; inactive; preserveAspect]
Now we try to do something more reasonable.  We use cddmathlink fuction AllVertices:
:[font = input; preserveAspect; startGroup]
?AllVertices
:[font = print; inactive; preserveAspect; endGroup]
AllVertices[m,d+1,A] generates all extreme points (vertices) and extreme
   rays of the convex polyhedron in R^(d+1) given as the solution set to an
   inequality system  A x >= 0 where  A is an m*(d+1) matrix  and 
   x=(1,x1,...,xd).  The output is {{extlist, linearity}, inclist} where
   extlist is  the extreme point list and inclist is the incidence list. 
   Each vertex (ray) has the first component 1 (0).  If the convex
   polyhedron is nonempty and has no vertices, extlist is a (nonunique) set
   of generators of the polyhedron where those generators in the linearity
   list are considered as linearity space (of points satisfying A (0, x1,
   x2, ...., xd) = 0)  generators.
:[font = subsection; inactive; preserveAspect; startGroup]
We can then compute ALL possibilities for cost at most, say One Dollar.
:[font = input; preserveAspect]
BudgetLimit=100;
:[font = input; preserveAspect; startGroup]
MatrixForm[dietproblem2=Append[dietproblem1, 
  {BudgetLimit, -3, -24, -13, -9, -20, -19}]]
:[font = output; output; inactive; preserveAspect]
MatrixForm[{{0, 1, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0}, 
 
   {0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 1, 0}, 
 
   {0, 0, 0, 0, 0, 0, 1}, {4, -1, 0, 0, 0, 0, 0}, {3, 0, -1, 0, 0, 0, 0}, 
 
   {2, 0, 0, -1, 0, 0, 0}, {8, 0, 0, 0, -1, 0, 0}, 
 
   {2, 0, 0, 0, 0, -1, 0}, {2, 0, 0, 0, 0, 0, -1}, 
 
   {-2000, 110, 205, 160, 160, 420, 260}, {-55, 4, 32, 13, 8, 4, 14}, 
 
   {-800, 2, 12, 54, 285, 22, 80}, {100, -3, -24, -13, -9, -20, -19}}]
;[o]
0       1       0       0       0       0       0

0       0       1       0       0       0       0

0       0       0       1       0       0       0

0       0       0       0       1       0       0

0       0       0       0       0       1       0

0       0       0       0       0       0       1

4       -1      0       0       0       0       0

3       0       -1      0       0       0       0

2       0       0       -1      0       0       0

8       0       0       0       -1      0       0

2       0       0       0       0       -1      0

2       0       0       0       0       0       -1

-2000   110     205     160     160     420     260

-55     4       32      13      8       4       14

-800    2       12      54      285     22      80

100     -3      -24     -13     -9      -20     -19
:[font = input; preserveAspect; startGroup]
{m2,d2}=Dimensions[dietproblem2]
:[font = output; output; inactive; preserveAspect; endGroup]
{16, 7}
;[o]
{16, 7}
:[font = input; preserveAspect; endGroup]
{{extlist,linearity},inclist}=AllVertices[m2,d2,Flatten[dietproblem2]];
:[font = input; preserveAspect; startGroup]
Length[extlist]
:[font = output; output; inactive; preserveAspect; endGroup]
17
;[o]
17
:[font = input; preserveAspect]
vlist=Map[Drop[#,1]&, extlist];
:[font = input; preserveAspect]
allsolutions=Union[Map[Prepend[#, N[c.#,3]]&, N[vlist,5]]];
:[font = input; preserveAspect; startGroup]
TableForm[table2=Prepend[allsolutions,var]]
:[font = output; output; inactive; preserveAspect; fontLeading = 0; endGroup; endGroup]
TableForm[{{"", "Oatmeal", "Chicken", "Eggs", "Milk", "Cherry Pie", "Pork\
 
     Beans"}, {92.5, 4., 0, 0, 4.5, 2., 0}, 
 
   {97.3333333333333, 4., 0, 0, 8., 0.6666666666666663, 0}, 
 
   {98.6035889070147, 4., 0, 0, 2.232952691680261, 2., 1.39510603588907}, 
 
   {100., 4., 0, 0, 8., 0.4172661870503596, 0.4028776978417263}, 
 
   {100., 4., 0, 0.4955752212389376, 8., 0.47787610619469, 0}, 
 
   {100., 4., 0.1872909698996654, 0, 8., 0.5752508361204011, 0}, 
 
   {100., 4., 0.6015037593984967, 0, 3.729323308270675, 2., 0}, 
 
   {100., 1.647058823529411, 0, 0, 6.117647058823531, 2., 0}, 
 
   {100., 2.808510638297873, 0, 0, 8., 0.978723404255319, 0}, 
 
   {100., 4., 0, 0, 2.179657768651608, 1.882819986310746, 
 
    1.617193702943191}, {100., 4., 0, 0, 2.209158679446219, 2., 
 
    1.479872204472844}, {100., 4., 0, 0, 5.333333333333334, 2., 0}, 
 
   {100., 4., 0, 0, 8., 0.7999999999999997, 0}, 
 
   {100., 4., 0, 1.025149700598807, 2.212215568862275, 2., 
 
    0.777005988023949}, {100., 4., 0, 1.875000000000002, 
 
    2.624999999999998, 2., 0}, 
 
   {100., 4., 0.1655305777133195, 0, 2.268813149625332, 2., 
 
    1.242523567802755}, {100., 3.741506869998489, 0, 0, 
 
    2.198037143288539, 2., 1.525955005284615}}]
;[o]
       Oatmeal   Chicken   Eggs      Milk     Cherry Pie   Pork Beans

92.5   4.        0         0         4.5      2.           0

97.3   4.        0         0         8.       0.66667      0

98.6   4.        0         0         2.233    2.           1.3951

100.   4.        0         0         8.       0.41727      0.40288

100.   4.        0         0.49558   8.       0.47788      0

100.   4.        0.18729   0         8.       0.57525      0

100.   4.        0.6015    0         3.7293   2.           0

100.   1.6471    0         0         6.1176   2.           0

100.   2.8085    0         0         8.       0.97872      0

100.   4.        0         0         2.1797   1.8828       1.6172

100.   4.        0         0         2.2092   2.           1.4799

100.   4.        0         0         5.3333   2.           0

100.   4.        0         0         8.       0.8          0

100.   4.        0         1.0251    2.2122   2.           0.77701

100.   4.        0         1.875     2.625    2.           0

100.   4.        0.16553   0         2.2688   2.           1.2425

100.   3.7415    0         0         2.198    2.           1.526
:[font = text; inactive; preserveAspect; endGroup]
The list is complete in the sense that any feasible menu of cost at most One Dollar is a combination of these seventeen (extreme) solutions.  One can find menus with Chicken, Eggs or Pork that might be much more desireble than the optimal menu.   Also it shows you cannot avoid Oatmeal nor Cherry pie within this budget to satisfy the nutritional needs.
:[font = section; inactive; Cclosed; preserveAspect; startGroup]
Disconnecting  cddmathlink
:[font = input; preserveAspect; startGroup]
Uninstall["~/Math/cddmathlink"]
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
"~/Math/cddmathlink"
;[o]
~/Math/cddmathlink
^*)
